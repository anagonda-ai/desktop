from Bio import SeqIO
import numpy as np
import networkx as nx

# Define distance metric function (replace with your chosen metric)
def distance_metric(cover_precentage, max_cover_precentage):
  # In this example, lower score indicates higher distance (less similarity)
  return cover_precentage / max_cover_precentage  # Normalize score between 0 and 1

# Load BLAST results
all_blast_input = input("Enter path to blast file:\n")


# Initialize distance matrix (assuming gene identifiers from FASTA headers)
ids_file = input("Enter path to ids file:\n")
gene_ids = [record.id for record in SeqIO.parse(ids_file,"fasta")]
max_cover_precentage = 0  # Initialize max_cover_precentage

# Find maximum bit score across all alignments and HSPs    
with open(all_blast_input, 'r') as blast_file:
  for line in blast_file:
    fields = line.strip().split("\t")
    current_cover_precentage = float(fields[2])
    max_cover_precentage = max(max_cover_precentage, current_cover_precentage)
    
# Now that max_cover_precentage is found, proceed with filling the distance matrix
distance_matrix = np.zeros((len(gene_ids), len(gene_ids)))  # Initialize with infinity

# Fill the matrix based on BLAST results (using the global max_cover_precentage)
with open(all_blast_input, 'r') as blast_file:
  for line in blast_file:
    fields = line.strip().split("\t")
    query_id, hit_id, cover_precentage = fields[0], fields[1], float(fields[2])
    distance = distance_metric(cover_precentage, max_cover_precentage)
    distance_matrix[gene_ids.index(query_id), gene_ids.index(hit_id)] = distance

# Create network from distance matrix
G = nx.DiGraph(distance_matrix)

# Calculate degree centrality
degree_centrality = nx.degree_centrality(G)

# Identify hub genes based on threshold or percentile
threshold = np.percentile(list(degree_centrality.values()), 99)
hub_genes = [f"{gene_ids[gene]}:{centrality}" for gene, centrality in degree_centrality.items() if centrality >= threshold]

with open("my_data.txt", "w") as file:
  file.write(f"Potential hub genes (degree centrality):\n")
  for hub in hub_genes:
    file.write(f"{hub}\n")
# print(f"Potential hub genes (degree centrality): {hub_genes}")