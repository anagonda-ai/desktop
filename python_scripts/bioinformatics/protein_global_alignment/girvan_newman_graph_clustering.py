import pandas as pd
import pickle
import networkx as nx
import numpy as np
import psutil
from sklearn.preprocessing import MinMaxScaler
from concurrent.futures import ThreadPoolExecutor

def load_and_normalize_scores(csv_file, threshold, chunk_size=5000):
    """Load similarity scores from CSV in chunks and normalize them using MinMaxScaler."""
    edges = []
    scores = []

    # Pass 1: Collect scores for normalization
    print("Scanning file for min/max scores...")
    for chunk in pd.read_csv(csv_file, chunksize=chunk_size):
        scores.extend(chunk["Similarity"].astype(float).tolist())

    # Fit MinMaxScaler
    scaler = MinMaxScaler()
    scores = np.array(scores).reshape(-1, 1)
    scaled_scores = scaler.fit_transform(scores).flatten()

    # Pass 2: Load edges with normalized scores
    print("Loading and normalizing scores using MinMaxScaler...")
    index = 0
    for chunk in pd.read_csv(csv_file, chunksize=chunk_size):
        chunk["Similarity"] = scaled_scores[index:index + len(chunk)]
        chunk = chunk[chunk["Similarity"] > threshold]  # Apply threshold
        edges.extend(list(zip(chunk["File1"], chunk["File2"], chunk["Similarity"])))
        index += len(chunk)

    return edges

# ----------- Step 2: Build Graph in Parallel ----------- #
def build_filtered_graph(edges):
    """Build a NetworkX graph from edges efficiently."""
    print("Building graph...")
    G = nx.Graph()
    
    with ThreadPoolExecutor() as executor:
        executor.map(lambda edge: G.add_edge(edge[0], edge[1], weight=edge[2]), edges)

    return G

# ----------- Step 3: Apply Clustering ----------- #
def apply_clustering(G):
    """Cluster the graph using Girvan-Newman algorithm in parallel."""
    if G.number_of_edges() == 0:
        return []

    communities_generator = nx.algorithms.community.girvan_newman(G)
    clusters = [list(cluster) for cluster in next(communities_generator)]  # First partition
    return clusters

# ----------- Step 4: Find Cluster Centers ----------- #
def find_cluster_centers(G, clusters):
    """Find the most connected node in each cluster."""
    centers = {}
    for i, cluster in enumerate(clusters):
        subgraph = G.subgraph(cluster)
        center = max(subgraph.degree, key=lambda x: x[1])[0]  # Highest-degree node
        centers[f"Cluster {i+1}"] = center
    return centers

# ----------- Step 5: Monitor Memory Usage ----------- #
def check_memory():
    """Check available system memory."""
    mem = psutil.virtual_memory()
    print(f"Memory Usage: {mem.percent}% ({mem.used / (1024 ** 3):.2f} GB used, {mem.available / (1024 ** 3):.2f} GB available)")
    return mem.percent < 90  # Stop if memory exceeds 90%

# ----------- Main Function ----------- #
def main(csv_file, threshold, output_graph, output_clusters, chunk_size=5_000_000):
    print("Starting pipeline...")
    
    if not check_memory():
        print("⚠️ Memory usage too high! Aborting.")
        return
    
    edges = load_and_normalize_scores(csv_file, threshold, chunk_size)

    if not check_memory():
        print("⚠️ Memory usage too high after loading scores! Aborting.")
        return
    
    G = build_filtered_graph(edges)
    print(f"Graph built with {len(G.nodes)} nodes and {len(G.edges)} edges.")

    if not check_memory():
        print("⚠️ Memory usage too high after graph construction! Aborting.")
        return
    
    print("Applying clustering...")
    clusters = apply_clustering(G)
    print(f"Clusters found: {len(clusters)}")

    print("Finding cluster centers...")
    cluster_centers = find_cluster_centers(G, clusters)

    # Save the graph
    with open(output_graph, 'wb') as f:
        pickle.dump(G, f)

    # Save clusters
    with open(output_clusters, 'w') as f:
        for i, cluster in enumerate(clusters):
            f.write(f"Cluster {i+1}: {', '.join(cluster)}\n")
        f.write("\nCluster Centers:\n")
        for cluster, center in cluster_centers.items():
            f.write(f"{cluster} center: {center}\n")

    print("Processing complete. Clusters and centers saved.")

if __name__ == "__main__":
    csv_file = "/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/protein_global_alignment/graph_output.csv"
    threshold = 0.5  # Adjust as needed
    output_graph = "/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/protein_global_alignment/girvan_newman_clustering_output/filtered_graph.pkl"
    output_clusters = "/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/protein_global_alignment/girvan_newman_clustering_output/clusters.txt"

    main(csv_file, threshold, output_graph, output_clusters)
