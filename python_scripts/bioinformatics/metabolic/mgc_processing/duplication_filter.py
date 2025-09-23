import os
import csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import networkx as nx
import subprocess
import pandas as pd
import time
import concurrent.futures
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
# Constants for BLAST parameters
EVALUE_THRESHOLD = 0.001
IDENTITY_THRESHOLD = 70.0
COVERAGE_THRESHOLD = 70.0

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

def parse_blast_results(blast_file_path):
    """
    Parse BLAST results and filter by identity threshold.
    Returns list of valid hits for graph construction.
    """
    hits = []
    
    try:
        with open(blast_file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue
                
                qseqid, sseqid, pident, length, evalue, bitscore = parts[:6]
                
                try:
                    pident = float(pident)
                    evalue = float(evalue)
                    length = int(length)
                    bitscore = float(bitscore)
                except ValueError:
                    continue
                
                # Filter by identity threshold
                if pident >= IDENTITY_THRESHOLD:
                    hits.append({
                        "query": qseqid,
                        "subject": sseqid,
                        "identity": pident,
                        "length": length,
                        "evalue": evalue,
                        "bitscore": bitscore
                    })
    
    except FileNotFoundError:
        print(f"⚠️  BLAST file not found: {blast_file_path}")
        return []
    
    return hits

def solve_maximum_weight_orthogonal_partitioning(query_genes, blast_hits):
    """
    Solve Maximum Weight Bipartite Matching using the Hungarian Algorithm.
    
    This is the OPTIMAL solution, not greedy.
    Finds the assignment of query genes to subject genes that maximizes total weight.
    """
    from collections import defaultdict
    
    # Try to import the best available solver
    solver = None
    try:
        from scipy.optimize import linear_sum_assignment
        import numpy as np
        solver = "scipy"
    except ImportError:
        try:
            import networkx as nx
            solver = "networkx"
        except ImportError:
            print("   ❌ Error: Need either scipy or networkx for optimal matching")
            return None, None
    
    # Build adjacency lists for bipartite graph
    gene_to_subjects = defaultdict(list)
    
    for hit in blast_hits:
        query_gene = hit["query"]
        if query_gene in query_genes:
            gene_to_subjects[query_gene].append({
                "subject": hit["subject"],
                "score": hit["bitscore"],
                "identity": hit["identity"],
                "evalue": hit["evalue"]
            })
    
    # Verify all query genes have potential matches
    for query_gene in query_genes:
        if query_gene not in gene_to_subjects or len(gene_to_subjects[query_gene]) == 0:
            print(f"   ❌ Query gene '{query_gene}' has no valid matches (≥{IDENTITY_THRESHOLD}% identity)")
            return None, None
    
    # Get all unique subject genes that appear in valid matches
    all_subjects = set()
    for subjects_list in gene_to_subjects.values():
        for subject_info in subjects_list:
            all_subjects.add(subject_info["subject"])
    
    if len(all_subjects) < len(query_genes):
        print(f"   ❌ Insufficient subjects: {len(all_subjects)} < {len(query_genes)} queries")
        return None, None
    
    query_list = sorted(list(query_genes))  # Sort for reproducibility
    subject_list = sorted(list(all_subjects))
    
    print(f"   🔥 Solving OPTIMAL Maximum Weight Bipartite Matching")
    print(f"      Queries: {len(query_list)}, Subjects: {len(subject_list)}")
    print(f"      Algorithm: {solver.upper()}")
    
    if solver == "scipy":
        return _hungarian_scipy(query_list, subject_list, gene_to_subjects)
    else:
        return _max_weight_matching_networkx(query_list, subject_list, gene_to_subjects)


def _hungarian_scipy(query_list, subject_list, gene_to_subjects):
    """
    Hungarian Algorithm implementation using scipy.
    Guaranteed to find the optimal maximum weight assignment.
    """
    import numpy as np
    from scipy.optimize import linear_sum_assignment
    
    n_queries = len(query_list)
    n_subjects = len(subject_list)
    
    # Create cost matrix
    # We use NEGATIVE weights because scipy minimizes, but we want to maximize
    cost_matrix = np.full((n_queries, n_subjects), np.inf)  # inf = impossible match
    
    # Fill matrix with actual costs (negative bitscores)
    for i, query_gene in enumerate(query_list):
        if query_gene in gene_to_subjects:
            for subject_info in gene_to_subjects[query_gene]:
                subject_gene = subject_info["subject"]
                j = subject_list.index(subject_gene)
                cost_matrix[i, j] = -subject_info["score"]  # Negative for maximization
    
    # Apply Hungarian algorithm
    print(f"      🧮 Running Hungarian algorithm on {n_queries}x{n_subjects} matrix...")
    
    try:
        row_indices, col_indices = linear_sum_assignment(cost_matrix)
    except Exception as e:
        print(f"   ❌ Hungarian algorithm failed: {e}")
        return None, None
    
    # Extract solution
    total_score = 0.0
    final_matching = {}
    
    print(f"      📋 Optimal assignment:")
    
    for i, j in zip(row_indices, col_indices):
        cost = cost_matrix[i, j]
        
        # Check if this assignment is valid (not infinite cost)
        if np.isinf(cost):
            print(f"   ❌ No valid complete matching exists")
            return None, None
        
        query_gene = query_list[i]
        subject_gene = subject_list[j]
        actual_score = -cost  # Convert back to positive
        
        # Find the original hit data
        hit_data = None
        for subject_info in gene_to_subjects[query_gene]:
            if subject_info["subject"] == subject_gene:
                hit_data = subject_info
                break
        
        if hit_data is None:
            print(f"   ❌ Internal error: lost hit data for {query_gene} → {subject_gene}")
            return None, None
        
        final_matching[query_gene] = hit_data
        total_score += actual_score
        
        print(f"         {query_gene} → {subject_gene} (score: {actual_score:.1f})")
    
    # Verify completeness
    if len(final_matching) != len(query_list):
        print(f"   ❌ Incomplete assignment: {len(final_matching)}/{len(query_list)}")
        return None, None
    
    print(f"   ✅ OPTIMAL solution found! Total score: {total_score:.2f}")
    return total_score, final_matching


def _max_weight_matching_networkx(query_list, subject_list, gene_to_subjects):
    """
    Maximum Weight Bipartite Matching using NetworkX.
    Alternative optimal algorithm when scipy is not available.
    """
    import networkx as nx
    
    # Create bipartite graph
    B = nx.Graph()
    
    # Add nodes with bipartite attribute (required for bipartite matching)
    query_nodes = {f"Q_{gene}": gene for gene in query_list}
    subject_nodes = {f"S_{gene}": gene for gene in subject_list}
    
    B.add_nodes_from(query_nodes.keys(), bipartite=0)  # Left side
    B.add_nodes_from(subject_nodes.keys(), bipartite=1)  # Right side
    
    # Add weighted edges for valid matches
    edge_count = 0
    for query_gene in query_list:
        query_node = f"Q_{query_gene}"
        for subject_info in gene_to_subjects[query_gene]:
            subject_gene = subject_info["subject"]
            subject_node = f"S_{subject_gene}"
            
            if subject_node in subject_nodes:
                B.add_edge(query_node, subject_node, weight=subject_info["score"])
                edge_count += 1
    
    print(f"      🕸️  Created bipartite graph: {len(query_nodes)} + {len(subject_nodes)} nodes, {edge_count} edges")
    
    # Find maximum weight matching
    try:
        # maxcardinality=True ensures we try to match as many nodes as possible
        matching = nx.algorithms.matching.max_weight_matching(B, maxcardinality=True)
        print(f"      ⚖️  Found matching with {len(matching)} edges")
        
    except Exception as e:
        print(f"   ❌ NetworkX matching failed: {e}")
        return None, None
    
    # Verify we have a complete matching for all queries
    matched_query_nodes = set()
    for edge in matching:
        node1, node2 = edge
        if node1.startswith('Q_'):
            matched_query_nodes.add(node1)
        elif node2.startswith('Q_'):
            matched_query_nodes.add(node2)
    
    if len(matched_query_nodes) < len(query_nodes):
        print(f"   ❌ Incomplete matching: {len(matched_query_nodes)}/{len(query_nodes)} queries matched")
        return None, None
    
    # Extract the solution
    total_score = 0.0
    final_matching = {}
    
    print(f"      📋 Optimal assignment:")
    
    for edge in matching:
        node1, node2 = edge
        
        # Determine which is query and which is subject
        if node1.startswith('Q_'):
            query_node, subject_node = node1, node2
        else:
            query_node, subject_node = node2, node1
        
        query_gene = query_nodes[query_node]
        subject_gene = subject_nodes[subject_node]
        
        # Find the hit data for this match
        hit_data = None
        for subject_info in gene_to_subjects[query_gene]:
            if subject_info["subject"] == subject_gene:
                hit_data = subject_info
                break
        
        if hit_data is None:
            print(f"   ❌ Internal error: lost hit data for {query_gene} → {subject_gene}")
            return None, None
        
        final_matching[query_gene] = hit_data
        total_score += hit_data["score"]
        
        print(f"         {query_gene} → {subject_gene} (score: {hit_data['score']:.1f})")
    
    print(f"   ✅ OPTIMAL solution found! Total score: {total_score:.2f}")
    return total_score, final_matching


def calculate_cluster_homology_score(query_cluster_file, target_cluster_file, fasta_map, all_blast_results):
    """
    Calculate homology score between two clusters using Maximum Weight Orthogonal Partitioning.
    
    Returns:
    - (score, matching_details) if clusters are homologous (all query genes matched orthogonally)
    - (None, None) if not homologous
    """
    
    # Get genes from query cluster
    query_fasta = fasta_map[query_cluster_file]
    query_genes = set()
    
    try:
        with open(query_fasta, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    gene_id = line.strip()[1:].split(" | ")[1]
                    query_genes.add(gene_id)
    except:
        return None, None
    
    if not query_genes:
        return None, None
    
    # Find corresponding BLAST results file
    query_base = os.path.splitext(os.path.basename(query_fasta))[0]
    target_base = os.path.splitext(os.path.basename(fasta_map[target_cluster_file]))[0]
    blast_filename = f"{query_base}.fasta_vs_{target_base}.fasta.txt"
    
    # Get BLAST hits for this cluster pair
    cluster_blast_hits = []
    for result in all_blast_results:
        if result["file"] == blast_filename:
            cluster_blast_hits.append(result)
    
    if not cluster_blast_hits:
        print(f"   ⚠️  No BLAST results found for {blast_filename}")
        return None, None
    
    print(f"🔍 Analyzing: {os.path.basename(query_cluster_file)} → {os.path.basename(target_cluster_file)}")
    print(f"   Query genes: {len(query_genes)}, BLAST hits: {len(cluster_blast_hits)}")
    
    # Solve Maximum Weight Orthogonal Partitioning
    total_score, matching = solve_maximum_weight_orthogonal_partitioning(query_genes, cluster_blast_hits)
    
    if total_score is None:
        print(f"❌ NOT HOMOLOGOUS: {os.path.basename(query_cluster_file)} → {os.path.basename(target_cluster_file)}")
        return None, None
    
    print(f"✅ HOMOLOGOUS: {os.path.basename(query_cluster_file)} → {os.path.basename(target_cluster_file)}")
    print(f"   Orthogonal score: {total_score:.2f}")
    
    return total_score, matching


def build_homology_graph_with_mwop(fasta_map, all_blast_results, blast_output_dir):
    """
    Build homology graph using Maximum Weight Orthogonal Partitioning.
    Edge exists from A→B if ALL genes in A can be matched orthogonally to distinct genes in B.
    """
    
    # Compute cluster sizes (gene counts)
    cluster_sizes = {}
    print("📏 Computing cluster sizes...")
    
    for cluster_file, fasta_path in fasta_map.items():
        try:
            with open(fasta_path, 'r') as f:
                gene_count = sum(1 for line in f if line.startswith('>'))
                cluster_sizes[cluster_file] = gene_count
                print(f"   {os.path.basename(cluster_file)}: {gene_count} genes")
        except:
            cluster_sizes[cluster_file] = 0
    
    # Build directed homology graph
    homology_graph = nx.DiGraph()
    homology_graph.add_nodes_from(fasta_map.keys())
    
    print(f"\n🧬 Building homology graph with Maximum Weight Orthogonal Partitioning...")
    
    processed_pairs = set()
    edges_added = 0
    homology_details = []
    
    for query_cluster in fasta_map.keys():
        for target_cluster in fasta_map.keys():
            if query_cluster == target_cluster:
                continue
            
            # Avoid duplicate work
            pair_key = (query_cluster, target_cluster)
            if pair_key in processed_pairs:
                continue
            processed_pairs.add(pair_key)
            
            # Test homology: query → target
            score_q_to_t, matching_q_to_t = calculate_cluster_homology_score(
                query_cluster, target_cluster, fasta_map, all_blast_results
            )
            
            # Test homology: target → query
            score_t_to_q, matching_t_to_q = calculate_cluster_homology_score(
                target_cluster, query_cluster, fasta_map, all_blast_results
            )
            
            # Add edges based on homology direction and cluster sizes
            size_q = cluster_sizes.get(query_cluster, 0)
            size_t = cluster_sizes.get(target_cluster, 0)
            
            # Case 1: Both directions are homologous
            if score_q_to_t is not None and score_t_to_q is not None:
                # Add edges from smaller to larger
                if size_q <= size_t:
                    homology_graph.add_edge(query_cluster, target_cluster, 
                                          weight=score_q_to_t, 
                                          type="bidirectional_homology",
                                          matching_count=len(matching_q_to_t))
                    edges_added += 1
                
                if size_t <= size_q:
                    homology_graph.add_edge(target_cluster, query_cluster, 
                                          weight=score_t_to_q, 
                                          type="bidirectional_homology",
                                          matching_count=len(matching_t_to_q))
                    edges_added += 1
                
                print(f"🔄 Bidirectional homology: {os.path.basename(query_cluster)} ↔ {os.path.basename(target_cluster)}")
                
                homology_details.append({
                    "query": os.path.basename(query_cluster),
                    "target": os.path.basename(target_cluster),
                    "type": "bidirectional",
                    "score_q_to_t": score_q_to_t,
                    "score_t_to_q": score_t_to_q
                })
            
            # Case 2: Only one direction is homologous (smaller → larger)
            elif score_q_to_t is not None and size_q <= size_t:
                homology_graph.add_edge(query_cluster, target_cluster, 
                                      weight=score_q_to_t, 
                                      type="unidirectional_homology",
                                      matching_count=len(matching_q_to_t))
                edges_added += 1
                print(f"➡️  Unidirectional homology: {os.path.basename(query_cluster)} → {os.path.basename(target_cluster)}")
                
                homology_details.append({
                    "query": os.path.basename(query_cluster),
                    "target": os.path.basename(target_cluster),
                    "type": "unidirectional",
                    "score": score_q_to_t
                })
            
            elif score_t_to_q is not None and size_t <= size_q:
                homology_graph.add_edge(target_cluster, query_cluster, 
                                      weight=score_t_to_q, 
                                      type="unidirectional_homology",
                                      matching_count=len(matching_t_to_q))
                edges_added += 1
                print(f"➡️  Unidirectional homology: {os.path.basename(target_cluster)} → {os.path.basename(query_cluster)}")
                
                homology_details.append({
                    "query": os.path.basename(target_cluster),
                    "target": os.path.basename(query_cluster),
                    "type": "unidirectional",
                    "score": score_t_to_q
                })
    
    print(f"\n🔗 Homology graph complete: {edges_added} homology edges found")
    
    # Find clusters and identify centrums
    print("🎯 Identifying centrums in homology clusters...")
    
    unique_representatives = []
    cluster_info = []
    
    for component in nx.weakly_connected_components(homology_graph):
        if len(component) == 1:
            centrum = next(iter(component))
            unique_representatives.append(centrum)
            cluster_info.append({
                "centrum": centrum,
                "cluster_size": 1,
                "members": [centrum],
                "is_valid_centrum": True,
                "centrum_type": "singleton"
            })
        else:
            # Find centrum: node that receives edges from ALL other nodes in component
            subgraph = homology_graph.subgraph(component)
            component_list = list(component)
            expected_in_degree = len(component_list) - 1
            
            valid_centrum = None
            max_in_degree = -1
            fallback_centrum = None
            
            for node in component_list:
                in_degree = subgraph.in_degree(node)
                if in_degree > max_in_degree:
                    max_in_degree = in_degree
                    fallback_centrum = node
                
                # True centrum: receives edges from ALL other nodes
                if in_degree == expected_in_degree:
                    valid_centrum = node
                    break
            
            if valid_centrum:
                unique_representatives.append(valid_centrum)
                cluster_info.append({
                    "centrum": valid_centrum,
                    "cluster_size": len(component_list),
                    "members": component_list,
                    "is_valid_centrum": True,
                    "centrum_type": "true_centrum"
                })
                print(f"🎯 True centrum: {os.path.basename(valid_centrum)} (receives {expected_in_degree} edges)")
            else:
                unique_representatives.append(fallback_centrum)
                cluster_info.append({
                    "centrum": fallback_centrum,
                    "cluster_size": len(component_list),
                    "members": component_list,
                    "is_valid_centrum": False,
                    "centrum_type": "fallback_centrum",
                    "inbound_edges": max_in_degree
                })
                print(f"⚠️  Fallback centrum: {os.path.basename(fallback_centrum)} ({max_in_degree}/{expected_in_degree} edges)")
    
    # Save comprehensive results
    dedup_file = os.path.join(blast_output_dir, "deduplicated_file_list.txt")
    cluster_details_file = os.path.join(blast_output_dir, "mwop_cluster_analysis.txt")
    homology_details_file = os.path.join(blast_output_dir, "homology_relationships.csv")
    
    with open(dedup_file, 'w') as f:
        f.write('\n'.join(unique_representatives) + '\n')
    
    # Save homology details as CSV
    if homology_details:
        df_homology = pd.DataFrame(homology_details)
        df_homology.to_csv(homology_details_file, index=False)
    
    # Generate detailed analysis report
    valid_centrums = sum(1 for info in cluster_info if info["is_valid_centrum"])
    total_clusters = len(cluster_info)
    
    with open(cluster_details_file, 'w') as f:
        f.write(f"Maximum Weight Orthogonal Partitioning (MWOP) Analysis\n{'=' * 70}\n\n")
        f.write(f"Algorithm Details:\n")
        f.write(f"- Identity threshold: ≥{IDENTITY_THRESHOLD}%\n")
        f.write(f"- All query genes must be matched orthogonally\n")
        f.write(f"- One-to-one matching (no gene reuse)\n")
        f.write(f"- Maximum weight selection\n")
        f.write(f"- Directed edges from smaller to larger clusters\n\n")
        
        f.write(f"Results Summary:\n")
        f.write(f"- Total input clusters: {len(fasta_map)}\n")
        f.write(f"- Homology relationships found: {len(homology_details)}\n")
        f.write(f"- Homology graph edges: {edges_added}\n")
        f.write(f"- Final clusters after merging: {total_clusters}\n")
        f.write(f"- True centrums: {valid_centrums}\n")
        f.write(f"- Fallback centrums: {total_clusters - valid_centrums}\n")
        f.write(f"- Representative clusters: {len(unique_representatives)}\n\n")
        
        # Detailed cluster information
        for i, info in enumerate(cluster_info, 1):
            f.write(f"Cluster {i}:\n")
            f.write(f"  Centrum: {os.path.basename(info['centrum'])}\n")
            f.write(f"  Type: {info['centrum_type']}\n")
            f.write(f"  Size: {info['cluster_size']} original clusters\n")
            f.write(f"  Gene counts: {[cluster_sizes.get(member, 0) for member in info['members']]}\n")
            f.write(f"  Valid centrum: {info['is_valid_centrum']}\n")
            if not info['is_valid_centrum']:
                f.write(f"  Inbound edges: {info.get('inbound_edges', 'N/A')}/{len(info['members'])-1}\n")
            f.write(f"  Members: {[os.path.basename(m) for m in info['members']]}\n\n")
    
    print(f"\n🏁 Maximum Weight Orthogonal Partitioning Complete!")
    print(f"   📊 Original clusters: {len(fasta_map)}")
    print(f"   🔗 Homology edges: {edges_added}")
    print(f"   🎯 Final clusters: {total_clusters}")
    print(f"   ✅ Representatives: {len(unique_representatives)}")
    print(f"   📄 Results: {dedup_file}")
    print(f"   📋 Analysis: {cluster_details_file}")
    print(f"   📈 Homology data: {homology_details_file}")
    
    return unique_representatives


def parse_and_filter_blast(blast_path):
    """Parse individual BLAST file and return filtered results"""
    results = parse_blast_results(blast_path)
    # Add filename to each result
    for result in results:
        result["file"] = os.path.basename(blast_path)
    return results

def postprocess(blast_output_dir, merged_list, temp_dir):
    print("\n🧩 Running postprocessing with Maximum Weight Orthogonal Partitioning...")
    results_dir = os.path.join(blast_output_dir, "results_fixed")
    summary_path = os.path.join(blast_output_dir, "blast_summary.csv")
    
    # Check if BLAST summary already exists
    if os.path.exists(summary_path):
        print(f"📋 Loading existing BLAST summary: {summary_path}")
        try:
            df = pd.read_csv(summary_path)
            all_blast_results = df.to_dict('records')
            print(f"✅ Loaded {len(all_blast_results)} BLAST hits from cache")
        except Exception as e:
            print(f"⚠️  Failed to load existing summary ({e}), reprocessing...")
            all_blast_results = _parse_blast_files(results_dir, summary_path)
    else:
        print("📖 BLAST summary not found, parsing all files...")
        all_blast_results = _parse_blast_files(results_dir, summary_path)

    # Prepare FASTA mapping
    fasta_map = {}
    for file_path in merged_list:
        fasta = prepare_fasta(file_path, temp_dir)
        if fasta:
            fasta_map[file_path] = fasta

    # Apply Maximum Weight Orthogonal Partitioning
    build_homology_graph_with_mwop(fasta_map, all_blast_results, blast_output_dir)


def _parse_blast_files(results_dir, summary_path):
    """Parse BLAST result files and save summary."""
    all_blast_results = []
    
    print("📖 Parsing BLAST results...")
    file_count = 0
    
    for root, _, files in os.walk(results_dir):
        for file in files:
            if file.endswith(".txt"):
                path = os.path.join(root, file)
                filtered = parse_and_filter_blast(path)
                all_blast_results.extend(filtered)
                file_count += 1
                print(f"   Parsed {len(filtered)} hits from {file}")
    
    print(f"📊 Processed {file_count} BLAST files, found {len(all_blast_results)} total hits")
    
    # Save summary of all BLAST results
    if all_blast_results:
        df = pd.DataFrame(all_blast_results)
        df.to_csv(summary_path, index=False)
        print(f"✅ Saved BLAST summary: {summary_path}")
    else:
        print("⚠️  No BLAST hits found!")
    
    return all_blast_results

def blast_all_vs_all_slurm(merged_list, output_root):
    slurm_script = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/blast_with_params.sh"
    blast_output_dir = os.path.join(output_root, "remove_duplicates_blast")
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
        
        while True:
            squeue_cmd = ["squeue", "-u", "alongonda", "--partition", "itaym"]
            squeue_result = subprocess.run(squeue_cmd, capture_output=True, text=True)
            running_jobs = len(squeue_result.stdout.strip().splitlines())
            if running_jobs < 150:
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
        squeue_cmd = "squeue -u alongonda | grep blast_"
        squeue_result = subprocess.run(squeue_cmd, shell=True, capture_output=True, text=True)
        running_jobs = len(squeue_result.stdout.strip().splitlines()) if squeue_result.stdout.strip() else 0
        if running_jobs == 0:
            break
        print(f"⏳ Waiting: {running_jobs} blast jobs still running. Sleeping 30s...")
        time.sleep(30)

def main():
    candidate_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test"
    mgc_dir = "/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/fasta_files"
    merged_list_file = os.path.join(candidate_dir, "merged_list.txt")
    blast_output_dir = os.path.join(candidate_dir, "remove_duplicates_blast")
    temp_dir = os.path.join(blast_output_dir, "temp_fastas")
    
    if not os.path.exists(merged_list_file):
        candidate_files = [
            os.path.join(candidate_dir, f) for f in os.listdir(mgc_dir)
            if f.endswith(".fasta")
        ]
        mgc_files = [
            os.path.join(mgc_dir, f) for f in os.listdir(mgc_dir)
            if f.endswith(".fasta")
        ]
        merged_list = candidate_files + mgc_files

        with open(merged_list_file, 'w') as f:
            for path in merged_list:
                f.write(f"{path}\n")
    else:
        with open(merged_list_file, 'r') as f:
            merged_list = [line.strip() for line in f]

    # blast_all_vs_all_slurm(merged_list, candidate_dir)
    # wait_till_blast_finish()
    postprocess(blast_output_dir, merged_list, temp_dir)

if __name__ == "__main__":
    main()