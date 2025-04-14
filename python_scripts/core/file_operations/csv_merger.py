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
                    df['pathway'] = file_name.split('_')[0]  # Extract pathway name from file name
                    dataframes.append(df)
                    print(f"Read CSV file: {file_path}")

    
    # Concatenate all dataframes
    merged_df = pd.concat(dataframes, ignore_index=True)
    
    # Save the merged dataframe to the given path
    merged_df.to_csv(output_csv_path, index=False)
    print(f"Saved merged CSV to: {output_csv_path}")

def main():
    csv_dir = "/groups/itay_mayrose/alongonda/datasets/KEGG_fasta_updated_fixed"
    output_csv_path = os.path.join(csv_dir, "merged_pathways.csv")
    merge_csv_files(csv_dir, output_csv_path)

if __name__ == "__main__":
    main()
