import os
import csv
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import tempfile
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
        fasta_out = file_path.replace(".csv", ".fasta")
        csv_to_fasta(file_path, fasta_out)
        return fasta_out
    else:
        return None

def build_bipartite_graph_fixed(candidate_fasta, mgc_fasta, temp_dir):
    """
    Fixed version that properly handles sequence IDs and creates a bipartite graph
    """
    db_prefix = os.path.join(temp_dir, "blast_db")
    blast_output = os.path.join(temp_dir, "blast_results.txt")
    
    # Make BLAST database from MGC sequences
    subprocess.run(["makeblastdb", "-in", mgc_fasta, "-dbtype", "prot", "-out", db_prefix], 
                   capture_output=True, check=True)
    
    # Run BLAST
    subprocess.run([
        "blastp", "-query", candidate_fasta, "-db", db_prefix,
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
        "-evalue", str(EVALUE_THRESHOLD), "-num_threads", str(NUM_THREADS), 
        "-out", blast_output
    ], capture_output=True, check=True)
    
    # Build bipartite graph from BLAST results
    graph = nx.Graph()
    
    # Add all nodes first
    for record in SeqIO.parse(candidate_fasta, "fasta"):
        graph.add_node(f"candidate_{record.id}", bipartite=0)
    for record in SeqIO.parse(mgc_fasta, "fasta"):
        graph.add_node(f"mgc_{record.id}", bipartite=1)
    
    # Add edges based on BLAST hits
    with open(blast_output) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                qid, sid, pident, length, evalue, bitscore = parts[:6]
                pident = float(pident)
                evalue = float(evalue)
                bitscore = float(bitscore)
                
                if pident >= IDENTITY_THRESHOLD and evalue <= EVALUE_THRESHOLD:
                    graph.add_edge(f"candidate_{qid}", f"mgc_{sid}", weight=bitscore)
    
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
            all_match_found = False

            for csv_filename in os.listdir(mgc_directory):
                if not csv_filename.endswith(".csv"):
                    continue
                
                csv_path = os.path.join(mgc_directory, csv_filename)
                
                # Convert CSV to temporary FASTA
                with tempfile.TemporaryDirectory() as temp_dir:
                    mgc_fasta = os.path.join(temp_dir, "mgc_temp.fasta")
                    csv_to_fasta(csv_path, mgc_fasta)
                    
                    # Build bipartite graph
                    graph = build_bipartite_graph_fixed(fasta_path, mgc_fasta, temp_dir)
                    
                    # Find maximum matching
                    matching = max_weight_matching(graph, maxcardinality=True)
                    
                    # Count candidate nodes to check if all matched
                    candidate_nodes = {n for n in graph.nodes() if n.startswith("candidate_")}
                    matched_candidates = {n for edge in matching for n in edge if n.startswith("candidate_")}
                    
                    if len(matched_candidates) == len(candidate_nodes) and len(candidate_nodes) > 0:
                        print(f"✅ All sequences in {fasta_path} matched to {csv_path}.")
                        out_file.write(f"The candidate {fasta_path} is matched to {csv_path} cluster.\n")
                        out_file.flush()
                        all_match_found = True
                        break

            if not all_match_found:
                print(f"❌ No complete matches found for {fasta_path} in any MGC CSV file.")
                out_file.write(f"Not all sequences in {fasta_path} have matches in any MGC CSV file.\n")
                out_file.flush()
                unmatched_candidates.append(fasta_path)
                
    return unmatched_candidates

def parse_and_filter_blast(blast_path):
    results = []
    with open(blast_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            qseqid, sseqid, pident, length, evalue, bitscore = parts[:6]
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
    """
    Efficient directed graph clustering implementation.
    Time complexity: O(E + V log V) where E = edges, V = vertices
    """
    
    # Pre-compute all file sizes once
    file_keys = {os.path.basename(v): k for k, v in fasta_map.items()}
    file_sizes = {}
    
    print("📏 Computing file sizes...")
    for file_path, fasta_path in fasta_map.items():
        try:
            # Fast line counting for FASTA files
            with open(fasta_path, 'r') as f:
                file_sizes[file_keys[os.path.basename(fasta_path)]] = sum(1 for line in f if line.startswith('>'))
        except:
            file_sizes[file_keys[os.path.basename(fasta_path)]] = 0
    
    # Single pass: build best scores lookup and directed graph simultaneously
    file_pair_best = defaultdict(lambda: {"score": -float('inf'), "data": None})
    file_graph = nx.DiGraph()
    
    # Add all nodes upfront
    file_graph.add_nodes_from(file_keys.values())
    
    print("🔍 Processing BLAST results...")
    for row in results:
        blast_file = row["file"]
        # Fix parsing: handle different filename formats
        if "_vs_" in blast_file:
            q_base, _, t_base = blast_file.partition("_vs_")
            t_base = t_base.replace(".txt", "")
        else:
            # Fallback parsing for different formats
            parts = blast_file.replace(".txt", "").split("_")
            if len(parts) >= 2:
                q_base = "_".join(parts[:-1])
                t_base = parts[-1]
            else:
                continue
        
        # Try different extensions
        q_file = None
        t_file = None
        for ext in [".fasta", ".fa"]:
            q_candidate = file_keys.get(q_base + ext)
            t_candidate = file_keys.get(t_base + ext)
            if q_candidate:
                q_file = q_candidate
            if t_candidate:
                t_file = t_candidate
        
        if q_file and t_file and q_file != t_file:
            # Keep only the best score between each pair
            pair_key = (q_file, t_file) if q_file < t_file else (t_file, q_file)
            
            if row["bitscore"] > file_pair_best[pair_key]["score"]:
                file_pair_best[pair_key] = {
                    "score": row["bitscore"],
                    "data": row
                }
    
    # Add directed edges efficiently
    print("🔗 Building directed graph...")
    edges_added = 0
    for (file1, file2), best_match in file_pair_best.items():
        size1 = file_sizes.get(file1, 0)
        size2 = file_sizes.get(file2, 0)
        
        # Direction: smaller -> larger (or lexicographic if equal)
        if size1 < size2 or (size1 == size2 and file1 < file2):
            source, target = file1, file2
        else:
            source, target = file2, file1
        
        file_graph.add_edge(source, target, 
                           weight=best_match["score"],
                           identity=best_match["data"]["identity"],
                           evalue=best_match["data"]["evalue"])
        edges_added += 1
    
    print(f"🔗 Added {edges_added} edges to graph")
    
    # Efficient clustering using weakly connected components
    print("🎯 Finding clusters and centrums...")
    unique_representatives = []
    cluster_info = []
    
    # Get weakly connected components (treats directed graph as undirected for connectivity)
    for component in nx.weakly_connected_components(file_graph):
        if len(component) == 1:
            centrum = next(iter(component))
            unique_representatives.append(centrum)
            cluster_info.append({
                "centrum": centrum,
                "cluster_size": 1,
                "members": [centrum],
                "is_valid_centrum": True
            })
        else:
            # Fast centrum detection using in-degree counting
            subgraph = file_graph.subgraph(component)
            component_list = list(component)
            expected_in_degree = len(component_list) - 1
            
            # Find node with in-degree equal to (cluster_size - 1)
            valid_centrum = None
            max_in_degree = -1
            fallback_centrum = None
            
            for node in component_list:
                in_degree = subgraph.in_degree(node)
                if in_degree > max_in_degree:
                    max_in_degree = in_degree
                    fallback_centrum = node
                
                if in_degree == expected_in_degree:
                    valid_centrum = node
                    break
            
            if valid_centrum:
                unique_representatives.append(valid_centrum)
                cluster_info.append({
                    "centrum": valid_centrum,
                    "cluster_size": len(component_list),
                    "members": component_list,
                    "is_valid_centrum": True
                })
            else:
                # Use fallback
                unique_representatives.append(fallback_centrum)
                cluster_info.append({
                    "centrum": fallback_centrum,
                    "cluster_size": len(component_list),
                    "members": component_list,
                    "is_valid_centrum": False,
                    "inbound_edges": max_in_degree
                })
    
    # Efficient file writing
    print("💾 Saving results...")
    dedup_file = os.path.join(blast_output_dir, "deduplicated_file_list.txt")
    cluster_details_file = os.path.join(blast_output_dir, "cluster_details.txt")
    
    # Write files efficiently
    with open(dedup_file, 'w') as f:
        f.write('\n'.join(unique_representatives) + '\n')
    
    # Generate summary statistics
    valid_centrums = sum(1 for info in cluster_info if info["is_valid_centrum"])
    total_clusters = len(cluster_info)
    
    with open(cluster_details_file, 'w') as f:
        # Write header
        f.write(f"Cluster Analysis Results\n{'=' * 50}\n\n")
        f.write(f"Total clusters: {total_clusters}\n")
        f.write(f"Valid centrums: {valid_centrums}\n")
        f.write(f"Invalid centrums: {total_clusters - valid_centrums}\n\n")
        
        # Write cluster details
        for i, info in enumerate(cluster_info, 1):
            f.write(f"Cluster {i}:\n")
            f.write(f"  Centrum: {info['centrum']}\n")
            f.write(f"  Size: {info['cluster_size']}\n")
            f.write(f"  Valid centrum: {info['is_valid_centrum']}\n")
            if not info['is_valid_centrum']:
                f.write(f"  Inbound edges: {info.get('inbound_edges', 'N/A')}\n")
            f.write(f"  Members: {', '.join(info['members'])}\n\n")
    
    print(f"\n🧬 Directed graph clustering complete:")
    print(f"   Original files: {len(file_keys)}")
    print(f"   Clusters found: {total_clusters}")
    print(f"   Valid centrums: {valid_centrums}")
    print(f"   Representatives: {len(unique_representatives)}")
    print(f"📄 Results saved to: {dedup_file}")
    print(f"📊 Cluster details saved to: {cluster_details_file}")
    
    return unique_representatives

def build_bipartite_graph_from_files(candidate_fasta, mgc_fasta):
    """
    Build bipartite graph using actual sequence IDs from FASTA files
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        db_prefix = os.path.join(temp_dir, "blast_db")
        blast_output = os.path.join(temp_dir, "blast_results.txt")
        
        # Make BLAST database
        subprocess.run(["makeblastdb", "-in", mgc_fasta, "-dbtype", "prot", "-out", db_prefix], 
                       capture_output=True, check=True)
        
        # Run BLAST
        subprocess.run([
            "blastp", "-query", candidate_fasta, "-db", db_prefix,
            "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
            "-evalue", str(EVALUE_THRESHOLD), "-num_threads", str(NUM_THREADS), 
            "-out", blast_output
        ], capture_output=True, check=True)
        
        # Build bipartite graph
        graph = nx.Graph()
        
        # Add candidate nodes
        for record in SeqIO.parse(candidate_fasta, "fasta"):
            graph.add_node(record.id, bipartite=0)
            
        # Add MGC nodes  
        for record in SeqIO.parse(mgc_fasta, "fasta"):
            graph.add_node(record.id, bipartite=1)
        
        # Add edges from BLAST results
        with open(blast_output) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    qid, sid, pident, length, evalue, bitscore = parts[:6]
                    pident = float(pident)
                    evalue = float(evalue)
                    bitscore = float(bitscore)
                    
                    if pident >= IDENTITY_THRESHOLD and evalue <= EVALUE_THRESHOLD:
                        graph.add_edge(qid, sid, weight=bitscore)
        
        return graph

def check_identity_fixed(mgc_directory, candidate_directory, output_file):
    """
    Fixed version of check_identity that properly builds bipartite graphs
    """
    unmatched_candidates = []
    with open(output_file, 'w') as out_file:
        for fasta_filename in os.listdir(candidate_directory):
            if not (fasta_filename.endswith(".fasta") or fasta_filename.endswith(".fa")):
                continue
            if "merged_mgc_candidartes.fasta" in fasta_filename:
                continue

            fasta_path = os.path.join(candidate_directory, fasta_filename)
            all_match_found = False

            # Count sequences in candidate file
            candidate_count = sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))

            for csv_filename in os.listdir(mgc_directory):
                if not csv_filename.endswith(".csv"):
                    continue
                
                csv_path = os.path.join(mgc_directory, csv_filename)
                
                # Convert CSV to temporary FASTA
                with tempfile.TemporaryDirectory() as temp_dir:
                    mgc_fasta = os.path.join(temp_dir, "mgc_temp.fasta")
                    csv_to_fasta(csv_path, mgc_fasta)
                    
                    # Build bipartite graph using actual sequence IDs
                    graph = build_bipartite_graph_from_files(fasta_path, mgc_fasta)
                    
                    # Find maximum matching
                    matching = max_weight_matching(graph, maxcardinality=True)
                    
                    # Count matched candidate sequences
                    matched_candidates = {n for edge in matching for n in edge 
                                        if graph.nodes[n].get('bipartite') == 0}
                    
                    print(f"🔍 {fasta_filename} vs {csv_filename}: {len(matched_candidates)}/{candidate_count} candidates matched")
                    
                    if len(matched_candidates) == candidate_count and candidate_count > 0:
                        print(f"✅ All sequences in {fasta_path} matched to {csv_path}.")
                        out_file.write(f"The candidate {fasta_path} is matched to {csv_path} cluster.\n")
                        out_file.flush()
                        all_match_found = True
                        break

            if not all_match_found:
                print(f"❌ No complete matches found for {fasta_path} in any MGC CSV file.")
                out_file.write(f"Not all sequences in {fasta_path} have matches in any MGC CSV file.\n")
                out_file.flush()
                unmatched_candidates.append(fasta_path)
                
    return unmatched_candidates

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
    
    # Use the fixed clustering function
    group_and_filter_redundant_files(results, fasta_map, blast_output_dir)

def blast_all_vs_all_slurm(merged_list, output_root):
    slurm_script = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/blast_with_params.sh"
    blast_output_dir = os.path.join(output_root, "blast_all_vs_all_fixed")
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

    for q_path, q_fasta in fasta_map.items():
        q_name = os.path.basename(q_fasta)
        q_id = os.path.splitext(q_name)[0].split('_')[-1] if '_' in q_name else os.path.splitext(q_name)[0].split('BGC')[-1]
        bucket = int(q_id) // 100
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
         
def wait_till_blast_finish():
    while True:
        # Fixed command - need to use shell=True for pipe operations
        squeue_cmd = "squeue -u alongonda | grep blast_"
        squeue_result = subprocess.run(squeue_cmd, shell=True, capture_output=True, text=True)
        running_jobs = len(squeue_result.stdout.strip().splitlines()) if squeue_result.stdout.strip() else 0
        if running_jobs == 0:
            break
        print(f"⏳ Waiting: {running_jobs} blast jobs still running. Sleeping 30s...")
        time.sleep(30)

def main():
    candidate_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test"
    mgc_dir = "/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/csv_files"
    merged_list_file = os.path.join(candidate_dir, "merged_list.txt")
    blast_output_dir = os.path.join(candidate_dir, "blast_all_vs_all_fixed")
    temp_dir = os.path.join(blast_output_dir, "temp_fastas")

    if not os.path.exists(merged_list_file):
        output_file = os.path.join(candidate_dir, "comparison_results.txt")
        # Use the fixed check_identity function
        unmatched_candidates = check_identity_fixed(mgc_dir, candidate_dir, output_file)

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
    wait_till_blast_finish()
    postprocess(blast_output_dir, merged_list, temp_dir)

if __name__ == "__main__":
    main()