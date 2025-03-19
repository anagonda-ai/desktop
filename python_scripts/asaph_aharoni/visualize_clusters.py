import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
import re

# Define file paths
MAPPING_FILE = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/output/dataset_organism_mapping.csv"
TREE_FILE = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/output/ncbi_species_tree_named.nwk"
COMPARISON_FILE = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/output/comparison_results.csv"

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
    
    def calculate_largest_chromosome_length(row):
        if not pd.isna(row['Largest Chromosome File']):
            largest_chromosome_file = row['Directory'] + '/' + row['Largest Chromosome File']
            if os.path.exists(largest_chromosome_file):
                chrom_df = pd.read_csv(largest_chromosome_file)
                chrom_df = chrom_df.sort_values(by='start')
                chrom_df["index"] = chrom_df["subject_gene"].apply(lambda x: x.split("|")[4] if isinstance(x, str) and len(x.split("|")) > 4 else None)
                length = chrom_df['end'].iloc[-1] - chrom_df['start'].iloc[0]
                index_diff = int(chrom_df['index'].iloc[-1]) - int(chrom_df['index'].iloc[0])
                return length, index_diff
        return None

    comparison_df[['Largest Chromosome Length', 'Index Difference']] = pd.DataFrame(
        comparison_df.apply(calculate_largest_chromosome_length, axis=1).tolist()
    )
    
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
    chromosome_lines = dict(zip(comparison_df['Organism'].replace(" ", "_"), comparison_df['Largest Chromosome Lines']))
    chromosome_cluster_length = dict(zip(comparison_df['Organism'].replace(" ", "_"), comparison_df['Largest Chromosome Length']))
    chromosome_cluster_length_genes = dict(zip(comparison_df['Organism'].replace(" ", "_"), comparison_df['Index Difference']))
    print("Sample Gene Counts:", list(gene_counts.items())[:5])
    print("Sample Largest Chromosome Lines:", list(chromosome_lines.items())[:5])
    print("Sample Largest Chromosome Length:", list(chromosome_cluster_length.items())[:5])
else:
    gene_counts = {}
    chromosome_lines = {}
    chromosome_cluster_length = {}
    chromosome_cluster_length_genes = {}
    
def label_function(clade):
    label = clade.name
    replaced_label = label.replace("_", " ") if label else label
    if label and replaced_label in gene_counts:
        gene_count = gene_counts.get(replaced_label, 0)
        chrom_lines = chromosome_lines.get(replaced_label, 0)
        dict_length = chromosome_cluster_length.get(replaced_label, 0)
        dict_length_genes = chromosome_cluster_length_genes.get(replaced_label, 0)
        chrom_length = dict_length if not pd.isna(dict_length) else 0
        chrom_length_genes = dict_length_genes if not pd.isna(dict_length_genes) else 0
        label = f"{label}:{gene_count}_GeneseMultiChromosome_{chrom_lines}_GeneseMaxChromosome_{chrom_length}bp_Length_{chrom_length_genes}_Genes"
    print(f"Labeling Clade: {label}")  # Debugging print
    return label

# Visualization
fig, ax = plt.subplots(figsize=(15, 15))
Phylo.draw(tree, axes=ax, label_func=label_function, do_show=False)

# Save the annotated tree

# **Step 6: Save tree in Newick format**
output_dir = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/output/"
os.makedirs(output_dir, exist_ok=True)
fig_path = os.path.join(output_dir, "phylogenetic_tree_with_cross_chromosome_genes.png")

plt.savefig(fig_path, dpi=300, bbox_inches="tight")
plt.show()

print("\n✅ Phylogenetic tree saved as 'phylogenetic_tree_with_cross_chromosome_genes.png'")
