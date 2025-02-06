import os
from Bio import SeqIO
from Bio import pairwise2
import networkx as nx
import concurrent.futures
import pickle

# Step 1: Parse each FASTA file to extract genes (sequences)
def parse_fasta(fasta_file):
    gene_sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_sequences.append(str(record.seq))
    return gene_sequences

# Step 2: Compute Pairwise Similarity using Needleman-Wunsch (Global Alignment)
def compute_similarity(fasta_file1, fasta_file2, genes1, genes2):
    similarity_scores = []

    # Perform global pairwise alignment for each gene in file1 with each gene in file2
    for gene1 in genes1:
        for gene2 in genes2:
            alignments = pairwise2.align.globalxx(gene1, gene2)
            best_alignment_score = alignments[0][2]  # Score is in the third column
            similarity_scores.append(best_alignment_score)

    avg_similarity_score = sum(similarity_scores) / len(similarity_scores) if similarity_scores else 0
    return fasta_file1, fasta_file2, avg_similarity_score

def process_results(future, G, output_path):
    try:
        file1, file2, similarity_score = future.result()
        G.add_edge(file1, file2, weight=similarity_score)
        print(f"Processed {file1} and {file2} with similarity score: {similarity_score:.4f}")
        print(f"Graph has {len(G.nodes)} nodes and {len(G.edges)} edges.")

        if len(G.edges) % 1000 == 0:
            nx.write_edgelist(G, output_path, delimiter=",", data=["weight"])
            print(f"Graph saved with {len(G.nodes)} nodes and {len(G.edges)} edges.")

    except Exception as e:
        print(f"Error processing a future: {e}")

# Step 3: Build the graph based on pairwise similarity scores using multiprocessing
def build_graph(fasta_files, output_path):
    G = nx.Graph()
    num_files = len(fasta_files)
    parsed_genes = {}

    # Parse all FASTA files once and cache results
    for fasta_file in fasta_files:
        parsed_genes[fasta_file] = parse_fasta(fasta_file)

    # Using ProcessPoolExecutor for parallel processing of file pairs
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = set()  # Using a set for better performance with removals

        for i in range(num_files):
            for j in range(i + 1, num_files):
                file1 = fasta_files[i]
                file2 = fasta_files[j]

                # Submit a task immediately
                future = executor.submit(compute_similarity, file1, file2, parsed_genes[file1], parsed_genes[file2])
                futures.add(future)

                # Process results as they complete
                done, _ = concurrent.futures.wait(futures, timeout=0, return_when=concurrent.futures.FIRST_COMPLETED)
                for completed_future in done:
                    process_results(completed_future, G, output_path)
                    futures.remove(completed_future)

    print(f"Graph built with {len(G.nodes)} nodes and {len(G.edges)} edges.")
    return G

# Step 4: Apply MWOP or a similar clustering algorithm
def apply_clustering(G):
    communities = nx.algorithms.community.girvan_newman(G)
    clusters = [list(cluster) for cluster in communities]
    print(f"Found {len(clusters)} clusters.")
    return clusters

# Main function to execute the workflow
def main(fasta_dir, output_path):
    # Step 1: Get all FASTA files in the directory
    fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if f.endswith(".fasta")]
    print(f"Found {len(fasta_files)} FASTA files in directory: {fasta_dir}")

    # Step 2: Build the graph from similarity scores
    print("Building the similarity graph...")
    G = build_graph(fasta_files, output_path)
    
    with open(output_path, 'wb') as f:
        pickle.dump(G, f)
        print(f"Graph saved with {len(G.nodes)} nodes and {len(G.edges)} edges.")

    # Step 3: Apply clustering (MWOP or modularity-based)
    print("Clustering files...")
    clusters = apply_clustering(G)

    # Output the clusters
    print("Clusters:")
    for idx, cluster in enumerate(clusters):
        print(f"Cluster {idx + 1}: {', '.join(cluster)}")


if __name__ == "__main__":
    # Provide the path to the directory containing FASTA files
    fasta_dir = "/groups/itay_mayrose/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/mgc_candidates_fasta_files_without_e2p2_filtered_test"  # Replace with your actual directory path
    output_path = "/groups/itay_mayrose/alongonda/desktop/graph_output.pkl"  # Specify the path where the graph will be saved incrementally

    # Execute the main function
    main(fasta_dir, output_path)
