import pandas as pd

# Load the CSV files
pathway_gene_count = pd.read_csv("/groups/itay_mayrose_nosnap/alongonda/plantcyc/all_organisms/pathway_gene_count.csv")
candidate_gene_count = pd.read_csv("/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/candidate_gene_counts.csv")

# Merge the DataFrames on the 'Pathway' column
merged_df = candidate_gene_count.merge(
    pathway_gene_count,
    on="Pathway",
    how="left",
    suffixes=("_candidate", "_pathway")
)

# Fill NaN values in case there are unmatched rows
merged_df["Gene Count_pathway"] = merged_df["Gene Count_pathway"].fillna(0).astype(int)

# Calculate the percentage of Gene Count_candidate in Gene Count_pathway
# Avoid division by zero and handle cases where percentage > 100
def calculate_percentage(row):
    if row["Gene Count_pathway"] == 0:
        return 0
    percentage = (row["Gene Count_candidate"] / row["Gene Count_pathway"]) * 100
    return "Error" if percentage > 100 else percentage

merged_df["Percentage"] = merged_df.apply(calculate_percentage, axis=1)

# Save the result to a new CSV or inspect
merged_df.to_csv("/groups/itay_mayrose_nosnap/alongonda/plantcyc/pmn_mgc_potential/mgc_candidates_process/candidate_gene_counts_with_pathway.csv", index=False)
print(merged_df)
