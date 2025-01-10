import pandas as pd

# Load the dataframes
df1 = pd.read_csv('/groups/itay_mayrose/alongonda/Plant_MGC/unique_clusters_sliding_window_outputs/unique_potential_groups_with_window_size_10.csv')
df2 = pd.read_csv('/groups/itay_mayrose/alongonda/Plant_MGC/unique_clusters_sliding_window_outputs_chromosome_sorted/unique_potential_groups_w10.csv')

# Find common rows
common_rows = pd.merge(df1, df2, how='inner')

# Find unique rows for each dataframe
unique_df1 = df1[~df1.apply(tuple,1).isin(df2.apply(tuple,1))]
unique_df2 = df2[~df2.apply(tuple,1).isin(df1.apply(tuple,1))]

# Save the results to CSV files
common_rows.to_csv('/groups/itay_mayrose/alongonda/desktop/python_scripts/pmn_scripts/sliding_windows/common_rows.csv', index=False)
unique_df1.to_csv('/groups/itay_mayrose/alongonda/desktop/python_scripts/pmn_scripts/sliding_windows/unique_rows_df1.csv', index=False)
unique_df2.to_csv('/groups/itay_mayrose/alongonda/desktop/python_scripts/pmn_scripts/sliding_windows/unique_rows_df2.csv', index=False)

# Display the results
print("Common Rows:")
print(common_rows)

print("\nUnique Rows in df1:")
print(unique_df1)

print("\nUnique Rows in df2:")
print(unique_df2)