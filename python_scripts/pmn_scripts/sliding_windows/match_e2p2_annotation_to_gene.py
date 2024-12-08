import pandas as pd

# Define the paths to the CSV files
csv_files = [
    '/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/plaza/plaza_genes.csv',
    '/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/ensembl/ensembl_genes.csv',
    '/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/phytozome/phytozome_genes.csv'
]

# Define the path to the DeepEC_Result.txt file
deep_ec_result_file = '/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/pmn.1732804635/log_files/3digit_EC_prediction.txt'

# Read the DeepEC_Result.txt file into a DataFrame
deep_ec_df = pd.read_csv(deep_ec_result_file, sep='\t')

# Initialize an empty DataFrame to hold all the data
combined_df = pd.DataFrame()

# Read and combine all CSV files
for csv_file in csv_files:
    df = pd.read_csv(csv_file, header=None)
    df.columns = ['Query ID', 'Species', 'Pathway']  # Assign column names to the DataFrame
    combined_df = pd.concat([combined_df, df])

# Remove duplicates based on 'Query ID'
combined_df = combined_df.drop_duplicates(subset='Query ID')

# Merge with the DeepEC_Result DataFrame
updated_df = combined_df.merge(deep_ec_df, on='Query ID', how='left')

# Drop rows without 'Predicted EC number'
updated_df = updated_df.dropna(subset=['Predicted EC number'])
# Drop rows where 'Predicted EC number' is 'EC number not predicted'
updated_df = updated_df[updated_df['Predicted EC number'] != 'EC number not predicted']

# Save the combined and updated DataFrame to a new CSV file
new_file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/pmn_mgc_potential/combined_updated_genes.csv'
updated_df.to_csv(new_file_path, index=False)