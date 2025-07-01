import os
import csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import networkx as nx
import subprocess
from networkx.algorithms.matching import max_weight_matching
import pandas as pd
import time

# Constants for BLAST parameters
EVALUE_THRESHOLD = 0.001
IDENTITY_THRESHOLD = 70.0
COVERAGE_THRESHOLD = 70.0
NUM_THREADS = 26

def read_sequences_from_csv(file_path):
    sequences = set()
    with open(file_path, mode='r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            sequence = row['Translation'].replace('*', '')
            sequences.add(sequence)
    return sequences

def read_sequences_from_fasta(file_path):
    sequences = set()
    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq).replace('*', '')
        sequences.add(sequence)
    return sequences

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

def prepare_fasta(file_path, temp_dir):
    if file_path.endswith(".fasta") or file_path.endswith(".fa"):
        return file_path
    elif file_path.endswith(".csv"):
        fasta_out = os.path.join(temp_dir, os.path.basename(file_path).replace(".csv", ".fasta"))
        csv_to_fasta(file_path, fasta_out)
        return fasta_out
    else:
        return None

def build_bipartite_graph(candidate_sequences, mgc_sequences):
    graph = nx.Graph()
    for candidate_seq in candidate_sequences:
        for mgc_seq in mgc_sequences:
            if candidate_seq in mgc_seq or mgc_seq in candidate_seq:
                graph.add_edge(candidate_seq, mgc_seq, weight=1)
    return graph

def check_identity(mgc_directory, candidate_directory, output_file):
    unmatched_candidates = []
    with open(output_file, 'w') as out_file:
        for fasta_filename in os.listdir(candidate_directory):
            if not (fasta_filename.endswith(".fasta") or fasta_filename.endswith(".fa")):
                continue
            if "merged_mgc_candidartes.fasta" in fasta_filename:
                continue

            fasta_path = os.path.join(candidate_directory, fasta_filename)
            candidate_sequences = read_sequences_from_fasta(fasta_path)
            all_match_found = False

            for csv_filename in os.listdir(mgc_directory):
                if not csv_filename.endswith(".csv"):
                    continue
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

def group_and_filter_redundant_files(results, fasta_map, blast_output_dir):
    file_graph = nx.Graph()
    file_keys = {os.path.basename(v): k for k, v in fasta_map.items()}
    for f in file_keys.values():
        file_graph.add_node(f)

    for row in results:
        blast_file = row["file"]
        q_base, _, t_base = blast_file.partition("_vs_")
        q_file = file_keys.get(q_base + ".fasta") or file_keys.get(q_base + ".fa")
        t_file = file_keys.get(t_base.replace(".txt", ".fasta")) or file_keys.get(t_base.replace(".txt", ".fa"))
        if q_file and t_file:
            file_graph.add_edge(q_file, t_file)

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

def postprocess(blast_output_dir, merged_list, temp_dir):
    print("\n🧩 Running postprocessing...")
    results_dir = os.path.join(blast_output_dir, "results_fixed")
    results = []

    fasta_map = {}
    for file_path in merged_list:
        fasta = prepare_fasta(file_path, temp_dir)
        if fasta:
            fasta_map[file_path] = fasta

    for root, _, files in os.walk(results_dir):
        for file in files:
            if file.endswith(".txt"):
                path = os.path.join(root, file)
                filtered = parse_and_filter_blast(path)
                results.extend(filtered)

    df = pd.DataFrame(results)
    summary_path = os.path.join(blast_output_dir, "blast_summary.csv")
    df.to_csv(summary_path, index=False)
    print(f"✅ Saved summary to {summary_path}")
    group_and_filter_redundant_files(results, fasta_map, blast_output_dir)

def blast_all_vs_all_slurm(merged_list, output_root):
    slurm_script = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/blast_with_params.sh"
    blast_output_dir = os.path.join(output_root, "blast_all_vs_all")
    results_dir = os.path.join(blast_output_dir, "results_fixed")
    blast_db_dir = os.path.join(blast_output_dir, "blast_dbs")
    temp_dir = os.path.join(blast_output_dir, "temp_fastas")

    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(blast_db_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)

    fasta_map = {}
    for file_path in merged_list:
        fasta = prepare_fasta(file_path, temp_dir)
        if fasta:
            fasta_map[file_path] = fasta
    
    can_run = False
    for q_path, q_fasta in fasta_map.items():
        q_name = os.path.basename(q_fasta)
        q_id = os.path.splitext(q_name)[0].split('_')[-1] if '_' in q_name else os.path.splitext(q_name)[0].split('BGC')[-1]
        bucket = int(q_id) // 100
        if bucket >= 32:
            can_run = True
        if not can_run:
            continue
        bucket_dir = os.path.join(results_dir, f"bucket_{bucket}")
        os.makedirs(bucket_dir, exist_ok=True)

        t_list = [t for t_path, t in fasta_map.items() if t_path != q_path]
        target_list_path = os.path.join(bucket_dir, f"{q_name}_targets.txt")
        with open(target_list_path, 'w') as f:
            for t in t_list:
                f.write(f"{t}\n")

        job_name = f"blast_{q_id}"
        # Wait until there are fewer than 100 running jobs for user 'alongonda'

        while True:
            squeue_cmd = ["squeue", "-u", "alongonda"]
            squeue_result = subprocess.run(squeue_cmd, capture_output=True, text=True)
            running_jobs = len(squeue_result.stdout.strip().splitlines())
            if running_jobs < 100:
                break
            print(f"⏳ Waiting: {running_jobs} jobs running for 'alongonda'. Sleeping 30s...")
            time.sleep(30)
        sbatch_cmd = [
            "sbatch",
            "--job-name", job_name,
            "--output", os.path.join(bucket_dir, f"{q_name}.out"),
            "--error", os.path.join(bucket_dir, f"{q_name}.err"),
            slurm_script,
            q_fasta,
            target_list_path,
            bucket_dir,
            blast_db_dir
        ]

        result = subprocess.run(sbatch_cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✅ Submitted SLURM job for {q_name}: {result.stdout.strip()}")
        else:
            print(f"❌ Failed to submit job for {q_name}: {result.stderr.strip()}")

def main():
    candidate_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/kegg_metabolic_output_g3_slurm/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test"
    mgc_dir = "/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/csv_files"
    merged_list_file = os.path.join(candidate_dir, "merged_list.txt")
    blast_output_dir = os.path.join(candidate_dir, "blast_all_vs_all")
    temp_dir = os.path.join(blast_output_dir, "temp_fastas")

    if not os.path.exists(merged_list_file):
        output_file = os.path.join(candidate_dir, "comparison_results.txt")
        unmatched_candidates = check_identity(mgc_dir, candidate_dir, output_file)

        mgc_files = [
            os.path.join(mgc_dir, f) for f in os.listdir(mgc_dir)
            if f.endswith(".csv")
        ]
        merged_list = unmatched_candidates + mgc_files

        with open(merged_list_file, 'w') as f:
            for path in merged_list:
                f.write(f"{path}\n")
    else:
        with open(merged_list_file, 'r') as f:
            merged_list = [line.strip() for line in f]

    blast_all_vs_all_slurm(merged_list, candidate_dir)
    postprocess(blast_output_dir, merged_list, temp_dir)

if __name__ == "__main__":
    main()