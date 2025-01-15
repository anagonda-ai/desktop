import os
import pandas as pd

def merge_csv_files(csv_dir, output_csv_path):
    # List to hold dataframes
    dataframes = []
    
    # Iterate through all files in the directory
    
    for root, dirs, files in os.walk(csv_dir):
        for file_name in files:
            if file_name.endswith('.csv'):
                file_path = os.path.join(root, file_name)
                if os.path.isfile(file_path):
                    df = pd.read_csv(file_path)
                    # Ensure 'id' column is treated as string
                    df['gene_name'] = df['gene_name'].astype(str)
                    dataframes.append(df)
                    print(f"Read CSV file: {file_path}")

    
    # Concatenate all dataframes
    merged_df = pd.concat(dataframes, ignore_index=True)
    
    # Ensure 'id' is string type before grouping
    merged_df['gene_name'] = merged_df['gene_name'].astype(str)
    
    # Generate the suffix for each id
    merged_df['suffix'] = merged_df.groupby('gene_name').cumcount()
    
    # Create the updated id with suffix
    merged_df['gene_name'] = merged_df['gene_name'] + '.' + merged_df['suffix'].astype(str)
    
    # Drop the temporary suffix column
    merged_df.drop(columns=['suffix'], inplace=True)
    
    # Save the merged dataframe to the given path
    merged_df.to_csv(output_csv_path, index=False)
    print(f"Saved merged CSV to: {output_csv_path}")

def main():
    csv_dir = "/groups/itay_mayrose/alongonda/full_genomes/ensembl/organisms"
    output_csv_path = "/groups/itay_mayrose/alongonda/full_genomes/annotations/ensembl_annotations.csv"
    merge_csv_files(csv_dir, output_csv_path)

if __name__ == "__main__":
    main()
