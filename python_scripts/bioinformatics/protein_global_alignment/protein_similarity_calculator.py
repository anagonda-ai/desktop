import os
import csv
import pickle
import networkx as nx
import concurrent.futures
import parasail
from Bio import SeqIO
from scipy.optimize import linear_sum_assignment
import psutil  # To check system memory
import numpy as np

# -------------- Step 1: Parse FASTA files ---------------- #
def parse_fasta(fasta_file):
    """Extract sequences from a FASTA file."""
    return [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]

# -------------- Step 2: Compute Global Alignment ---------------- #
def compute_similarity(seq1, seq2):
    """Compute global alignment score using parasail."""
    alignment = parasail.nw_stats_striped_16(seq1, seq2, 10, 1, parasail.blosum62)
    return alignment.score

# -------------- Step 4: Graph Builder ---------------- #
def graph_builder(G, output_csv):
    """Reads similarity results from CSV and builds the graph incrementally."""
    with open(output_csv, 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header
        for row in reader:
            file1, file2, score = row[0], row[1], float(row[2])
            G.add_edge(file1, file2, weight=score)
            
def process_batch(batch, output_csv):
    """Process multiple pairs in one batch to reduce I/O overhead."""
    results = []
    for fasta_file1, fasta_file2, genes1, genes2 in batch:
        similarity_matrix = np.array([[compute_similarity(g1, g2) for g2 in genes2] for g1 in genes1])
        row_ind, col_ind = linear_sum_assignment(similarity_matrix, maximize=True)
        max_sum = similarity_matrix[row_ind, col_ind].sum()
        
        num_pairs = len(row_ind)
        avg_max_sum = max_sum / num_pairs if num_pairs else 0
        results.append([fasta_file1, fasta_file2, avg_max_sum])
        print(f"Processed: {fasta_file1} vs {fasta_file2} with similarity {avg_max_sum}")
        

    # Write batch results
    with open(output_csv, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(results)

# -------------- Step 5: Build Graph with Immediate Future Removal ---------------- #
def build_graph(fasta_files, output_csv, output_graph, batch_size=100):
    """Computes similarity for all file pairs and builds the graph."""
    num_files = len(fasta_files)
    parsed_genes = {f: parse_fasta(f) for f in fasta_files}

    # Clear previous CSV file and add header
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["File1", "File2", "Similarity"])

    # Use multiprocessing for parallel processing
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = set()

        # 🔥 **NEW: Loops through all file pairs** 🔥
        for i in range(num_files):
            for j in range(i + 1, num_files):  # Ensures every possible (i, j) is covered
                batch = [(fasta_files[i], fasta_files[j], parsed_genes[fasta_files[i]], parsed_genes[fasta_files[j]])]

                # Submit batch for processing
                future = executor.submit(process_batch, batch, output_csv)
                futures.add(future)

                # Manage completed futures to avoid memory overload
                if len(futures) >= batch_size:
                    done, _ = concurrent.futures.wait(futures, return_when=concurrent.futures.FIRST_COMPLETED)
                    for completed_future in done:
                        completed_future.result()  # Ensure errors are caught
                        futures.remove(completed_future)

        # Wait for any remaining futures to finish
        concurrent.futures.wait(futures)

    # Build graph from the CSV
    G = nx.Graph()
    graph_builder(G, output_csv)

    # Save graph
    with open(output_graph, 'wb') as f:
        pickle.dump(G, f)

    return G

# -------------- Step 6: Apply Clustering ---------------- #
def apply_clustering(G):
    """Apply clustering algorithm to the similarity graph."""
    communities = nx.algorithms.community.girvan_newman(G)
    clusters = [list(cluster) for cluster in next(communities)]
    return clusters

# -------------- Main Function ---------------- #
def main(fasta_dir, output_csv, output_graph):
    fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if f.endswith(".fasta")]
    print(f"Found {len(fasta_files)} FASTA files.")

    print("Building similarity graph...")
    G = build_graph(fasta_files, output_csv, output_graph)

    print(f"Graph saved with {len(G.nodes)} nodes and {len(G.edges)} edges.")

    print("Applying clustering...")
    clusters = apply_clustering(G)

    print(f"Clusters found: {len(clusters)}")
    for i, cluster in enumerate(clusters):
        print(f"Cluster {i+1}: {', '.join(cluster)}")

if __name__ == "__main__":
    fasta_dir = "/groups/itay_mayrose/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/mgc_candidates_fasta_files_without_e2p2_filtered_test"
    output_csv = "/groups/itay_mayrose/alongonda/desktop/graph_output.csv"
    output_graph = "/groups/itay_mayrose/alongonda/desktop/graph_output.pkl"
    
    main(fasta_dir, output_csv, output_graph)
