import os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

def create_dataframe(file_path):
    # Read the file into a DataFrame, skipping initial comment lines
    df = pd.read_csv(file_path, sep='\t', comment='#', header=None)
    return df

def process_file(file_path, output_directory):
    df = create_dataframe(file_path)
    
    # Filter rows where the value in the 3rd column is "gene"
    df_filtered = df[df.iloc[:, 2] == "gene"]
    
    
    # Create a new DataFrame with the specified columns
    df_transformed = pd.DataFrame({
        'id': df_filtered.iloc[:, 8].str.split(";").str[1].str.replace(".2.0.227", "").str.replace("Name=", ""),
        'start': df_filtered.iloc[:, 3],
        'end': df_filtered.iloc[:, 4]
    })
    
    # Save the transformed DataFrame as a new CSV file
    transformed_output_file_path = os.path.join(output_directory, f"{os.path.splitext(os.path.basename(file_path))[0]}_transformed.csv")
    df_transformed.to_csv(transformed_output_file_path, index=False)
    print(f"Transformed and saved {file_path} to {transformed_output_file_path}")
    
    return df_filtered

def process_directory(directory, output_directory):
    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)
    
    dataframes = []
    file_paths = []
    
    # Iterate through each file in the directory
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".gene.gff3"):
                file_path = os.path.join(root, file)
                file_paths.append(file_path)
    
    # Use ThreadPoolExecutor to process files concurrently
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_file, file_path, output_directory) for file_path in file_paths]
        for future in as_completed(futures):
            dataframes.append(future.result())
    
    return dataframes

def main():
    directory = "/groups/itay_mayrose_nosnap/alongonda/full_genomes/phytozome_test/Phytozome"
    output_directory = "/groups/itay_mayrose_nosnap/alongonda/full_genomes/phytozome_test/processed_annotations"
    dataframes = process_directory(directory, output_directory)
    for df in dataframes:
        print("DataFrame:")
        print(df.head())  # Print the first few rows of each DataFrame

if __name__ == "__main__":
    main()