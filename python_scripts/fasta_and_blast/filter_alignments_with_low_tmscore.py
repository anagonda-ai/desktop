import pandas as pd
# Load the TSV file into a DataFrame
df = pd.read_csv('/groups/itay_mayrose/alongonda/desktop/single_protein_blast/arabidopsis_foldseek_calculations/aln_tmscore_final_folder/aln_tmscore_final_preprocesses.tsv', sep='\t', header=None)

# Filter out rows where the first column value ends with ".cif"
pdb_df = df[~df[0].str.endswith('.cif')]

# Find the indices of rows with the highest values in column 2 for each unique value in column 0
max_indices = pdb_df.groupby(0)[2].idxmax()

# Select the rows with the highest values in column 2 for each unique value in column 0
highest_values_df = pdb_df.loc[max_indices]

treshold = 0.8

filtered_df = highest_values_df[highest_values_df.iloc[:, 2] >= treshold]
# Now you can perform operations on the DataFrame, such as filtering or analyzing data
print(f"We Found Plantcyc matches with foldseek score >= {treshold} for {int(filtered_df.shape[0])}/{int(highest_values_df.shape[0])} Arabidopsis genes")  # Display the first few rows of the DataFrame

filtered_df.to_csv(f'/groups/itay_mayrose/alongonda/desktop/single_protein_blast/arabidopsis_foldseek_calculations/aln_tmscore_final_folder/filtered_tmscore_final_with_{treshold}_treshold.tsv', sep='\t', index=False, header=False)
