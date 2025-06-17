import os
import csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import networkx as nx
import subprocess
from Bio.Blast.Applications import NcbiblastpCommandline
import concurrent.futures
from networkx import Graph
from networkx.algorithms.matching import max_weight_matching
import pandas as pd

# Constants for BLAST parameters
EVALUE_THRESHOLD = 0.001
IDENTITY_THRESHOLD = 70.0
COVERAGE_THRESHOLD = 70.0
NUM_THREADS = 26

# Function to read sequences from a CSV file
def read_sequences_from_csv(file_path):
    sequences = set()
    with open(file_path, mode='r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            sequence = row['Translation'].replace('*', '')  # Remove '*' characters
            sequences.add(sequence)
    return sequences

# Function to read sequences from a FASTA file
def read_sequences_from_fasta(file_path):
    sequences = set()
    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq).replace('*', '')  # Remove '*' characters
        sequences.add(sequence)
    return sequences

# Function to build a bipartite graph from candidate and MGC sequences
def build_bipartite_graph(candidate_sequences, mgc_sequences):
    graph = Graph()
    for candidate_seq in candidate_sequences:
        for mgc_seq in mgc_sequences:
            if candidate_seq in mgc_seq or mgc_seq in candidate_seq:
                graph.add_edge(candidate_seq, mgc_seq, weight=1)
    return graph

# Function to check identity between candidate sequences and MGC sequences
def check_identity(mgc_directory, candidate_directory, output_file):
    unmatched_candidates = []
    with open(output_file, 'w') as out_file:
        for fasta_filename in os.listdir(candidate_directory):
            if fasta_filename.endswith(".fasta") or fasta_filename.endswith(".fa"):
                if fasta_filename.endswith("merged_mgc_candidartes.fasta"):
                    continue
                fasta_path = os.path.join(candidate_directory, fasta_filename)
                candidate_sequences = read_sequences_from_fasta(fasta_path)
                all_match_found = False
                
                for csv_filename in os.listdir(mgc_directory):
                    if csv_filename.endswith(".csv"):
                        csv_path = os.path.join(mgc_directory, csv_filename)
                        mgc_sequences = read_sequences_from_csv(csv_path)
                        
                        graph = build_bipartite_graph(candidate_sequences, mgc_sequences)
                        matching = max_weight_matching(graph, maxcardinality=True)
                        
                        if len(matching) == len(candidate_sequences):
                            out_file.write(f"The candidate {fasta_path} is matched to {csv_path} cluster.\n")
                            all_match_found = True
                            break
                
                if not all_match_found:
                    out_file.write(f"Not all sequences in {fasta_path} have matches in any MGC CSV file.\n")
                    unmatched_candidates.append(fasta_path)
    return unmatched_candidates

# Function to convert CSV to FASTA format
def csv_to_fasta(csv_file, output_fasta):
    records = []
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if 'Translation' in row and row['Translation']:
                sequence = row['Translation'].replace('*', '')
                gene_id = row.get('Gene') or row.get('gene') or f"seq_{len(records)}"
                record = SeqRecord(Seq(sequence), id=gene_id, description="")
                records.append(record)
    if records:
        SeqIO.write(records, output_fasta, "fasta")

# Function to prepare FASTA files from CSV or existing FASTA files
def prepare_fasta(file_path, temp_dir):
    if file_path.endswith(".fasta") or file_path.endswith(".fa"):
        return file_path
    elif file_path.endswith(".csv"):
        fasta_out = os.path.join(temp_dir, os.path.basename(file_path).replace(".csv", ".fasta"))
        csv_to_fasta(file_path, fasta_out)
        return fasta_out
    else:
        return None

# Function to run a single BLAST search    
def run_single_blast(query_fasta, target_fasta, db_dir, output_dir):
    
    # Extract the ID from the FASTA filename    
    def extract_id(fasta_path):
        return os.path.splitext(os.path.basename(fasta_path))[0].split('_')[-1] if '_' in os.path.basename(fasta_path) else os.path.splitext(os.path.basename(fasta_path))[0].split('BGC')[-1]

    query_id = extract_id(query_fasta)
    bucket = int(query_id) // 100
    bucket_dir = os.path.join(output_dir, f"bucket_{bucket}")
    os.makedirs(bucket_dir, exist_ok=True)

    output_path = os.path.join(bucket_dir, f"{os.path.basename(query_fasta)}_vs_{os.path.basename(target_fasta)}.txt")
    db_name = os.path.join(db_dir, os.path.basename(target_fasta).replace(".fasta", ""))
    temp_output = f"{output_path}.tmp"

    if os.path.exists(output_path):
        print(f"BLAST output already exists: {output_path}")
        return output_path

    print(f"Running BLAST for {os.path.basename(query_fasta)} against {os.path.basename(target_fasta)}")

    if os.path.getsize(query_fasta) == 0 or os.path.getsize(target_fasta) == 0:
        print(f"⚠️ Skipping BLAST: One of the FASTA files is empty: {query_fasta}, {target_fasta}")
        return None

    if not os.path.exists(f"{db_name}.pin"):
        result = subprocess.run(["makeblastdb", "-in", target_fasta, "-dbtype", "prot", "-out", db_name], capture_output=True, text=True)
        if result.returncode != 0 or not os.path.exists(f"{db_name}.pin"):
            print(f"❌ makeblastdb failed for {target_fasta}: {result.stderr}")
            return None

    blastp_cline = NcbiblastpCommandline(
        query=query_fasta,
        db=db_name,
        evalue=EVALUE_THRESHOLD,
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        out=output_path
    )

    try:
        stdout, stderr = blastp_cline()
    except subprocess.CalledProcessError as e:
        print(f"❌ BLAST failed for {query_fasta} vs {target_fasta}")
        print(f"Command: {blastp_cline}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        if os.path.exists(output_path):
            os.remove(output_path)
        return None

    print(f"BLAST completed: {output_path}")
    return output_path

# Function to parse and filter BLAST results
def parse_and_filter_blast(blast_path):
    results = []
    with open(blast_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            qseqid, sseqid, pident, length, *_ , evalue, bitscore = parts[:12]
            pident = float(pident)
            evalue = float(evalue)
            length = int(length)
            bitscore = float(bitscore)
            if pident >= IDENTITY_THRESHOLD and evalue <= EVALUE_THRESHOLD:
                results.append({
                    "query": qseqid,
                    "subject": sseqid,
                    "identity": pident,
                    "length": length,
                    "evalue": evalue,
                    "bitscore": bitscore,
                    "file": os.path.basename(blast_path)
                })
    return results

# Function to group and filter redundant files based on homologous hits
def group_and_filter_redundant_files(results, fasta_map, blast_output_dir):
    file_graph = nx.Graph()

    # Initialize all files as nodes
    file_keys = {os.path.basename(v): k for k, v in fasta_map.items()}
    for f in file_keys.values():
        file_graph.add_node(f)

    # Add edges between files with any homologous hits
    for row in results:
        blast_file = row["file"]
        q_base, _, t_base = blast_file.partition("_vs_")
        q_file = file_keys.get(q_base + ".fasta") or file_keys.get(q_base + ".fa")
        t_file = file_keys.get(t_base.replace(".txt", ".fasta")) or file_keys.get(t_base.replace(".txt", ".fa"))
        if q_file and t_file:
            file_graph.add_edge(q_file, t_file)

    # Group homologous files
    unique_representatives = []
    for component in nx.connected_components(file_graph):
        representative = sorted(component)[0]
        unique_representatives.append(representative)

    dedup_file = os.path.join(blast_output_dir, "deduplicated_file_list.txt")
    with open(dedup_file, 'w') as f:
        for rep in unique_representatives:
            f.write(f"{rep}\n")

    print(f"\n🧬 Removed homologous redundancy: reduced {len(file_keys)} → {len(unique_representatives)} files.")
    print(f"📄 Deduplicated file list saved to: {dedup_file}")
    return unique_representatives

# Function to run all vs all BLAST in parallel
def blast_all_vs_all_parallel(merged_list, output_root):
    blast_output_dir = os.path.join(output_root, "blast_all_vs_all")
    blast_output_dir_results = os.path.join(blast_output_dir, "results_fixed")
    blast_db_dir = os.path.join(blast_output_dir, "blast_dbs")
    temp_dir = os.path.join(blast_output_dir, "temp_fastas")
    os.makedirs(blast_output_dir, exist_ok=True)
    os.makedirs(blast_output_dir_results, exist_ok=True)
    os.makedirs(blast_db_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)

    # Convert all to FASTA and store map
    fasta_map = {}
    for file_path in merged_list:
        fasta = prepare_fasta(file_path, temp_dir)
        if fasta:
            fasta_map[file_path] = fasta

    # Build all pairwise combinations
    tasks = []
    for q_path, q_fasta in fasta_map.items():
        for t_path, t_fasta in fasta_map.items():
            if q_path != t_path:
                tasks.append((q_fasta, t_fasta))

    # Run BLASTs concurrently
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_THREADS) as executor:
        future_to_blast = {
            executor.submit(run_single_blast, q, t, blast_db_dir, blast_output_dir_results): (q, t)
            for q, t in tasks
        }
        for future in concurrent.futures.as_completed(future_to_blast):
            blast_file = future.result()
            filtered = parse_and_filter_blast(blast_file)
            results.extend(filtered)

    # Save all filtered hits into one summary file
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(blast_output_dir, "blast_summary.csv"), index=False)
    dedup_list = group_and_filter_redundant_files(results, fasta_map, blast_output_dir)

# Write the deduplicated FASTA files to the output directory
def main():
    mgc_directory = "/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/csv_files"
    candidate_directory = "/groups/itay_mayrose/alongonda/Plant_MGC/kegg_metabolic_output_w10_g3/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test"
    output_file = os.path.join(candidate_directory, "comparison_results.txt")
    
    if not os.path.exists(os.path.join(candidate_directory, "merged_list.txt")):
        # Check the identity of the candidates
        unmatched_candidates = check_identity(mgc_directory, candidate_directory, output_file)
        
        # Get the MGC filepaths
        mgc_filepaths = os.listdir(mgc_directory)
        # Prepend the full path to MGC filepaths
        mgc_filepaths = [os.path.join(mgc_directory, filename) for filename in mgc_filepaths]
        
        # Merge the unmatched candidates with the MGC filepaths
        merged_list = unmatched_candidates + mgc_filepaths
        
        # Write the merged list to a file
        with open(os.path.join(candidate_directory, "merged_list.txt"), 'w') as merged_file:
            for item in merged_list:
                merged_file.write(f"{item}\n")
    else:
        with open(os.path.join(candidate_directory, "merged_list.txt"), 'r') as merged_file:
            merged_list = [line.strip() for line in merged_file.readlines()]
        print(f"⚠️ Merged list already exists: {os.path.join(candidate_directory, 'merged_list.txt')}")
            
    blast_all_vs_all_parallel(merged_list, candidate_directory)

# Function to run the script
if __name__ == "__main__":
    main()