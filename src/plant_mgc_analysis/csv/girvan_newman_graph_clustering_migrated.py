
class GirvanNewmanGraphClusteringProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: girvan_newman_graph_clustering.py."""
    
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
    nodes = set()
    scores = []

    # Pass 1: Collect scores for normalization
    # logger.info("Scanning file for min/max scores...")
    # for chunk in pd.read_csv(csv_file, chunksize=chunk_size):
    #     scores.extend(chunk["Similarity"].astype(float).tolist())

    # Fit MinMaxScaler
    # scaler = MinMaxScaler()
    # scores = np.array(scores).reshape(-1, 1)
    # scaled_scores = scaler.fit_transform(scores).flatten()

    # Pass 2: Load edges with normalized scores
    logger.info("Loading and normalizing scores using MinMaxScaler...")
    index = 0
    for chunk in pd.read_csv(csv_file, chunksize=chunk_size):
        # chunk["Similarity"] = scaled_scores[index:index + len(chunk)]
        chunk["Similarity"] = chunk["Similarity"].astype(float)
        nodes.update(chunk["File1"].tolist())
        nodes.update(chunk["File2"].tolist())
        chunk = chunk[chunk["Similarity"] > threshold]  # Apply threshold
        edges.extend(list(zip(chunk["File1"], chunk["File2"], chunk["Similarity"])))
        index += len(chunk)

    return edges, nodes

# ----------- Step 2: Build Graph in Parallel ----------- #
def build_filtered_graph(edges, nodes):
    """Build a NetworkX graph from edges efficiently."""
    logger.info("Building graph...")
    G = nx.Graph()
    
    G.add_nodes_from(nodes)
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
    logger.info(f"Memory Usage: {mem.percent}% ({mem.used / (1024 ** 3):.2f} GB used, {mem.available / (1024 ** 3):.2f} GB available)")
    return mem.percent < 90  # Stop if memory exceeds 90%

def filter_nodes_by_centers(G, clusters, centers):
    """Filter nodes, keeping only the center and isolated nodes not connected to any cluster."""
    filtered_nodes = set()

    # Add the cluster centers to the filtered nodes
    for cluster in clusters:
        center = centers[f"Cluster {clusters.index(cluster) + 1}"]
        filtered_nodes.add(center)

    # Find isolated nodes (nodes with no edges) and add them if they are not part of any cluster
    isolated_nodes = set(G.nodes) - filtered_nodes
    for node in isolated_nodes:
        if len(list(G.neighbors(node))) == 0:
            filtered_nodes.add(node)

    # Filter the graph to keep only the filtered nodes
    G_filtered = G.subgraph(filtered_nodes).copy()
    return G_filtered

# ----------- Step 6: Save the Centers and Isolated Nodes ----------- #
def save_filtered_graph(G, clusters, centers, output_graph):
    """Save only the centers and isolated nodes to the graph file."""
    G_filtered = filter_nodes_by_centers(G, clusters, centers)
    
    # Save all nodes to a CSV file
    all_nodes = list(G_filtered.nodes)
    nodes_df = pd.DataFrame(all_nodes, columns=["Node"])
    nodes_csv_file = output_graph.replace(".pkl", "_nodes.csv")
    nodes_df.to_csv(nodes_csv_file, index=False)
    # Save the filtered graph
    with open(output_graph, 'wb') as f:
        pickle.dump(G_filtered, f)

    logger.info(f"Filtered graph saved with {len(G_filtered.nodes)} nodes and {len(G_filtered.edges)} edges.")


# ----------- Main Function ----------- #
def main(csv_file, threshold, output_graph, output_clusters, chunk_size=5_000_000):
    logger.info("Starting pipeline...")
    
    if not check_memory():
        logger.info("⚠️ Memory usage too high! Aborting.")
        return
    
    edges, nodes = load_and_normalize_scores(csv_file, threshold, chunk_size)
    logger.info(f"Scores loaded from CSV with {len(nodes)} nodes.")

    if not check_memory():
        logger.info("⚠️ Memory usage too high after loading scores! Aborting.")
        return
    
    G = build_filtered_graph(edges, nodes)
    logger.info(f"Graph built with {len(G.nodes)} nodes and {len(G.edges)} edges.")

    if not check_memory():
        logger.info("⚠️ Memory usage too high after graph construction! Aborting.")
        return
    
    logger.info("Applying clustering...")
    clusters = apply_clustering(G)
    logger.info(f"Clusters found: {len(clusters)}")

    logger.info("Finding cluster centers...")
    cluster_centers = find_cluster_centers(G, clusters)
    
    # Save filtered graph (only centers and isolated nodes)
    save_filtered_graph(G, clusters, cluster_centers, output_graph)

    logger.info("Processing complete. Clusters and centers saved.")

if __name__ == "__main__":
    csv_file = "/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/protein_global_alignment/graph_output_mwp.csv"
    threshold = 150  # Adjust as needed
    output_graph = "/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/protein_global_alignment/girvan_newman_clustering_output/filtered_graph_updated.pkl"
    output_clusters = "/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/protein_global_alignment/girvan_newman_clustering_output/clusters_updated.txt"

    main(csv_file, threshold, output_graph, output_clusters)
