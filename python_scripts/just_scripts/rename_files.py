import os
import pandas as pd

def read_species_info(species_info_path):
    # Read the species information CSV file into a DataFrame
    df = pd.read_csv(species_info_path, delimiter='\t')
    keys = df.iloc[:, 0].values
    values = df.iloc[:, 1].values
    # Create a dictionary mapping species codes to common names
    species_map = pd.Series(values, index=keys).to_dict()
    return species_map

def rename_files(directory, species_map):
    # Iterate through each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            # Get the base name of the file without the extension
            base_name = os.path.splitext(filename)[0]
            # Check if the base name is in the species map
            if base_name in species_map:
                file_name = species_map[base_name].lower().replace(" ", "_")
                full_name = f"{file_name}.csv"
                old_file_path = os.path.join(directory, filename)
                new_file_path = os.path.join(directory, full_name)
                os.rename(old_file_path, new_file_path)
                print(f"Renamed {filename} to {full_name}")

def main():
    directory = "/groups/itay_mayrose/alongonda/full_genomes/plaza/processed_annotations"
    species_info_path = "/groups/itay_mayrose/alongonda/desktop/species_information.csv"
    
    # Read the species information
    species_map = read_species_info(species_info_path)
    
    # Rename the files in the directory
    rename_files(directory, species_map)

if __name__ == "__main__":
    main()