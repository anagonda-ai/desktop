import os
import pandas as pd

def create_dataframe(file_path):
    # Read the file into a DataFrame, skipping the initial comment lines
    df = pd.read_csv(file_path, sep='\t', comment='#', header=0)
    
    # Extract the columns that start with 'GENE-NAME-' and 'GENE-ID-'
    gene_name_cols = [col for col in df.columns if col.startswith('GENE-NAME')]
    gene_id_cols = [col for col in df.columns if col.startswith('GENE-ID')]
    
    # Rename the columns to match the desired format
    df.columns = ['UNIQUE-ID', 'NAME'] + gene_name_cols + gene_id_cols
    
    return df

def process_pathways_files(base_dir):
    dataframes = []
    for root, _, files in os.walk(base_dir):
        for file_name in files:
            if file_name == "pathways.col":
                file_path = os.path.join(root, file_name)
                df = create_dataframe(file_path)
                print(f"Created DataFrame for: {file_path}")
                dataframes.append(df)
                
    # Concatenate all DataFrames into one
    merged_df = pd.concat(dataframes, ignore_index=True)
    return merged_df


def main():
    base_dir = "/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms"
    output_csv_path = "/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms/merged_pathways.csv"
    merged_df = process_pathways_files(base_dir)
    print("Merged DataFrame:")
    print(merged_df.head())  # Print the first few rows of the merged DataFrame
    
    # Save the merged DataFrame to a CSV file
    merged_df.to_csv(output_csv_path, index=False)
    print(f"Saved merged DataFrame to: {output_csv_path}")

if __name__ == "__main__":
    main()