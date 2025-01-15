import pandas as pd

# Load CSV files
candidate_gene_counts = pd.read_csv("/groups/itay_mayrose/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/results/candidate_gene_counts_with_pathway.csv")
unique_clusters = pd.read_csv("/groups/itay_mayrose/alongonda/Plant_MGC/unique_clusters_sliding_window_outputs_chromosome_sorted/unique_potential_groups_w10.csv")

# Extract the pathway name and occurrence for matching
candidate_gene_counts['Pathway_Extracted'] = candidate_gene_counts['Pathway (Occurrence)']
unique_clusters['Pathway_Extracted'] = unique_clusters['pathway']

# Merge the DataFrames on the extracted Pathway name
merged_df = pd.merge(
    candidate_gene_counts,
    unique_clusters[['Pathway_Extracted', 'pathway_cluster_genes', 'start', 'end']],
    how='left',
    on='Pathway_Extracted'
)

# Drop the temporary extracted column
merged_df = merged_df.drop(columns=['Pathway_Extracted'])

# Convert the distance to kilobase pairs (kbp) and add the unit as a string
merged_df['Distance'] = ((merged_df['end'] - merged_df['start']) / 1000).astype(str) + " kbp"

# Save the result to a new CSV
merged_df.to_csv("/groups/itay_mayrose/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/results/candidate_gene_counts_with_pathway_start_end.csv", index=False)
