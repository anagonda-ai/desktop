
class ClusterCounterProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: cluster_counter.py."""
    
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

def analyze_clusters(fasta_file):
  """
  Analyzes a FASTA file and calculates the number of clusters and genes per cluster.

  Args:
      fasta_file (str): Path to the FASTA file containing gene sequences.

  Returns:
      dict: A dictionary where keys are cluster names and values are lists of genes within that cluster.
  """
  clusters = {}
  current_cluster = None

  with open(fasta_file, 'r') as f:
    for line in f:
      line = line.rstrip()  # Remove trailing newline character

      if line.startswith('>'):
        # New cluster identifier
        cluster_name = line[1:].split('|')[0]
        current_cluster = cluster_name
        clusters.setdefault(cluster_name, [])  # Create empty list if cluster doesn't exist
      else:
        # Gene sequence, add to current cluster
        if current_cluster:
          clusters[current_cluster].append(line)
        else:
          # No cluster defined yet, skip this gene
          pass

  return clusters

# Example usage
fasta_file = input("Enter fasta file path:\n")  # Replace with your FASTA file path
cluster_data = analyze_clusters(fasta_file)

ordered_clusters = sorted(cluster_data.items(), key=lambda x: len(x[1]), reverse=True)

# Calculate total number of genes and clusters
total_genes = sum(len(cluster[1]) for cluster in ordered_clusters)
num_clusters = len(ordered_clusters)

# Calculate average cluster size (if there are clusters)
average_size = total_genes / num_clusters if num_clusters > 0 else 0

output_file = os.path.join(input("Enter output dir path:\n"), "cluster_count.txt")
with open(output_file, 'w') as f:
    f.write(f"Number of clusters: {len(ordered_clusters)}\n")
    f.write(f"Average cluster size: {average_size:.2f}\n")  # Format average with 2 decimals
    for cluster_name, genes in ordered_clusters:
      f.write(f"Cluster: {cluster_name}, Number of genes: {len(genes)}\n")
