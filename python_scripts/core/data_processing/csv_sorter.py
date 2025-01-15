import os
import pandas as pd

def reorder_csv_by_start(input_csv_dir, output_csv_dir):
    # Create the output directory if it doesn't exist
    os.makedirs(output_csv_dir, exist_ok=True)
    
    # Iterate through all files in the directory
    for file_name in os.listdir(input_csv_dir):
        file_path = os.path.join(input_csv_dir, file_name)
        if os.path.isfile(file_path) and file_name.endswith('.csv'):
            # Read the CSV file into a DataFrame
            df = pd.read_csv(file_path)
            
            # Sort the DataFrame by the 'start' column
            df_sorted = df.sort_values(by=['start', 'end'])
            
            # Save the sorted DataFrame to the new directory
            output_file_path = os.path.join(output_csv_dir, file_name)
            df_sorted.to_csv(output_file_path, index=False)
            print(f"Saved sorted CSV to: {output_file_path}")

def main():
    input_csv_dir = "/groups/itay_mayrose/alongonda/full_genomes/ensembl/processed_annotations"
    output_csv_dir = "/groups/itay_mayrose/alongonda/full_genomes/ensembl/processed_annotations_sorted"
    reorder_csv_by_start(input_csv_dir, output_csv_dir)

if __name__ == "__main__":
    main()