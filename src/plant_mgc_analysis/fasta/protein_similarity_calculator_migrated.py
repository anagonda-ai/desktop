
class ProteinSimilarityCalculatorProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: protein_similarity_calculator.py."""
    
    def __init__(self, **kwargs):
        """Initialize processor."""
        super().__init__(**kwargs)
        self.settings = get_settings()
    
    def validate_input(self, data):
        """Validate input data."""
        pass  # Implement validation
    
    def process(self, data, **kwargs):
        """Process data."""
        # Original script logic here
        pass


from ..config.settings import get_settings
from ..core.base import BioinformaticsProcessor
from ..core.exceptions import ProcessingError, ValidationError
from ..utils.file_operations import file_manager
from loguru import logger
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
        if abs(len(genes1) - len(genes2)) > 1:
            logger.info(f"Skipping: {fasta_file1} vs {fasta_file2} due to gene count difference")
            num_pairs = 0
        else:
            # Create the bipartite graph
            G = nx.Graph()

            # Compute similarity scores and add edges
            for i, g1 in enumerate(genes1):
                for j, g2 in enumerate(genes2):
                    score = compute_similarity(g1, g2)
                    if score > 0:  # Keep only meaningful edges
                        G.add_edge(g1, g2, weight=score)
            # Solve MWOP using maximum weight matching
            matching = nx.algorithms.matching.max_weight_matching(G, maxcardinality=False, weight="weight")
            num_pairs = len(matching)
            # Compute max sum
            max_sum = sum(G[u][v]["weight"] for u, v in matching)
        
        avg_max_sum = max_sum / num_pairs if not num_pairs == 0 else 0
        results.append([fasta_file1, fasta_file2, avg_max_sum])
        logger.info(f"Processed: {fasta_file1} vs {fasta_file2} with similarity {avg_max_sum}")
        

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
    logger.info(f"Found {len(fasta_files)} FASTA files.")

    logger.info("Building similarity graph...")
    G = build_graph(fasta_files, output_csv, output_graph)

    logger.info(f"Graph saved with {len(G.nodes)} nodes and {len(G.edges)} edges.")

    logger.info("Applying clustering...")
    clusters = apply_clustering(G)

    logger.info(f"Clusters found: {len(clusters)}")
    for i, cluster in enumerate(clusters):
        logger.info(f"Cluster {i+1}: {', '.join(cluster)}")

if __name__ == "__main__":
    fasta_dir = "/groups/itay_mayrose/alongonda/datasets/plantcyc/pmn_mgc_potential/mgc_candidates_process/mgc_candidates_fasta_files_without_e2p2_filtered_test"
    output_csv = "/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/protein_global_alignment/test.csv"
    output_graph = "/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/protein_global_alignment/test.pkl"
    
    main(fasta_dir, output_csv, output_graph)
