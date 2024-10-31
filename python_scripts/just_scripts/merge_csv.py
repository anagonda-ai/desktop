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
            # Ensure 'id' column is treated as string
            df['id'] = df['id'].astype(str)
            dataframes.append(df)
            print(f"Read CSV file: {file_path}")
    
    # Concatenate all dataframes
    merged_df = pd.concat(dataframes, ignore_index=True)
    
    # Ensure 'id' is string type before grouping
    merged_df['id'] = merged_df['id'].astype(str)
    
    # Generate the suffix for each id
    merged_df['suffix'] = merged_df.groupby('id').cumcount()
    
    # Create the updated id with suffix
    merged_df['id'] = merged_df['id'] + '.' + merged_df['suffix'].astype(str)
    
    # Drop the temporary suffix column
    merged_df.drop(columns=['suffix'], inplace=True)
    
    # Save the merged dataframe to the given path
    merged_df.to_csv(output_csv_path, index=False)
    print(f"Saved merged CSV to: {output_csv_path}")

def main():
    csv_dir = "/groups/itay_mayrose_nosnap/alongonda/full_genomes/annotations/"
    output_csv_path = "/groups/itay_mayrose_nosnap/alongonda/full_genomes/annotations/merged_annotations.csv"
    merge_csv_files(csv_dir, output_csv_path)

if __name__ == "__main__":
    main()
