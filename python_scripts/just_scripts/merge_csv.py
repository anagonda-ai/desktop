import os
import pandas as pd

def merge_csv_files(csv_dir, output_csv_path):
    # List to hold dataframes
    dataframes = []
    
    # Iterate through all files in the directory
    for file_name in os.listdir(csv_dir):
        file_path = os.path.join(csv_dir, file_name)
        if os.path.isfile(file_path) and file_name.endswith('.csv'):
            df = pd.read_csv(file_path)
            dataframes.append(df)
            print(f"Read CSV file: {file_path}")
    
    # Concatenate all dataframes
    merged_df = pd.concat(dataframes, ignore_index=True)
    
    # Save the merged dataframe to the given path
    merged_df.to_csv(output_csv_path, index=False)
    print(f"Saved merged CSV to: {output_csv_path}")

def main():
    
    csv_dir = "/groups/itay_mayrose_nosnap/alongonda/full_genomes/phytozome_test/processed_annotations"
    output_csv_path = "/groups/itay_mayrose_nosnap/alongonda/full_genomes/phytozome_test/phytozome_annotations.csv"
    merge_csv_files(csv_dir, output_csv_path)

if __name__ == "__main__":
    main()