import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
import re

# Define file paths
MAPPING_FILE = "/groups/itay_mayrose/alongonda/desktop/dataset_organism_mapping.csv"
TREE_FILE = "/groups/itay_mayrose/alongonda/desktop/ncbi_species_tree_named.nwk"
COMPARISON_FILE = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/blast_results_chromosome_separated/best_hits_by_organism/comparison_results.csv"

# Load phylogenetic tree
tree = Phylo.read(TREE_FILE, "newick")

# Load dataset organism mapping file
if os.path.exists(MAPPING_FILE):
    mapping_df = pd.read_csv(MAPPING_FILE)
    mapping_dict = dict(zip(mapping_df['Original Filename'], mapping_df['Organism']))
else:
    mapping_df = None
    mapping_dict = {}
    print("Warning: Dataset organism mapping file not found.")

# Load comparison results file
if os.path.exists(COMPARISON_FILE):
    comparison_df = pd.read_csv(COMPARISON_FILE)
    
    # Extract only the relevant part of the directory to match mapping_dict
    def clean_directory_name(directory):
        return directory.split("/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/blast_results_chromosome_separated/best_hits_by_organism/")[1].split(".gene_transformed_filtered")[0]  # Extract the first meaningful part
    
    comparison_df['Cleaned_Name'] = comparison_df['Directory'].apply(clean_directory_name)
    
    # Match Cleaned_Name with any substring in mapping_dict keys
    def find_matching_organism(cleaned_name):
        for key in mapping_dict:
            if cleaned_name in key:
                return mapping_dict[key]
        return 'Unknown'
    
    comparison_df['Organism'] = comparison_df['Cleaned_Name'].apply(find_matching_organism)
else:
    comparison_df = None
    print("Warning: Comparison results file not found.")

# Create a dictionary for organism to gene count mapping
if comparison_df is not None:
    gene_counts = dict(zip(comparison_df['Organism'].replace(" ", "_"), comparison_df['Cross Chromosome Lines']))  # Adjust column name if necessary
else:
    gene_counts = {}
    
def label_function(clade):
    label = clade.name
    replaced_label = label.replace("_", " ") if label else label
    if label and replaced_label in gene_counts:
        label = f"{label}:{gene_counts.get(replaced_label, 0)}"
    print(f"Labeling Clade: {label}")  # Debugging print
    return label

# Visualization
fig, ax = plt.subplots(figsize=(15, 15))
Phylo.draw(tree, axes=ax, label_func=label_function, do_show=False)

# Save the annotated tree
plt.savefig("phylogenetic_tree_with_cross_chromosome_genes.png", dpi=300, bbox_inches="tight")
plt.show()

print("\n✅ Phylogenetic tree saved as 'phylogenetic_tree_with_cross_chromosome_genes.png'")
