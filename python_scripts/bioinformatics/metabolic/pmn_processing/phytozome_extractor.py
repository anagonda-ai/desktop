import os
import re
import time
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing


def find_file_for_organism(base_path, organism_name, mode="phytozome"):
    """
    Recursively search for the file matching the organism name.
    
    :param base_path: Base directory to start searching
    :param organism_name: Organism name to match
    :param mode: 'phytozome' or 'ensembl' to determine file format
    :return: Full path to the matching file or None
    """
    for root, dirs, files in os.walk(base_path):
        if mode == "phytozome":
            for filename in files:
                if filename.endswith('.protein.fa') and organism_name in filename:
                    print(f"Found file for {organism_name}: {os.path.join(root, filename)}")
                    return os.path.join(root, filename)
        elif mode == "ensembl":
            if organism_name in dirs:
                organism_dir = os.path.join(root, organism_name)
                csv_file = os.path.join(organism_dir, "extracted_data.csv")
                if os.path.exists(csv_file):
                    print(f"Found file for {organism_name}: {csv_file}")
                    return csv_file
        elif mode == "plaza":
            if organism_name in dirs:
                organism_dir = os.path.join(root, organism_name)
                print(f"Organism dir: {organism_dir}")
                for filename in os.listdir(organism_dir):
                    print(f"Checking file: {filename}")
                    if not filename.startswith("start"):
                        file_path = os.path.join(organism_dir, filename)
                        print(f"Found file for {organism_name}: {file_path}")
                        return file_path
    print(f"WARNING: No file found for organism: {organism_name}")
    return None


def extract_gene_protein_sequence(file_path, gene_name, mode="phytozome"):
    """
    Extract protein sequence for a specific gene, based on the file format.

    :param file_path: Path to the input file
    :param gene_name: Gene identifier to search for
    :param mode: 'phytozome' (FASTA format) or 'ensembl' (CSV format)
    :return: Protein sequence or None if not found
    """
    try:
        if mode == "phytozome" or mode == "plaza":
            new_gene_name = gene_name.split(';')[1].split('=')[1] if mode == "phytozome" else gene_name
            # Handle FASTA format for phytozome
            with open(file_path, 'r') as f:
                found = False
                sequence_lines = []
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):  # Header line
                        if found:  # If already found, stop
                            break
                        if new_gene_name in line:
                            found = True
                    elif found:
                        sequence_lines.append(line)
                if sequence_lines:
                    sequence = ''.join(sequence_lines)
                    print(f"Found sequence for {gene_name} ({mode}): {sequence[:50]}...")
                    return sequence
                print(f"No sequence found for {gene_name} ({mode})")
                return None

        elif mode == "ensembl":
            # Handle CSV format for ensembl
            fixed_gene_name = gene_name.lower().replace("_p", "_t") # Fix gene name for Ensembl
            df = pd.read_csv(file_path)
            gene_row = df[df['transcript'].apply(lambda x: fixed_gene_name.find(x.lower()) != -1)]
            if not gene_row.empty:
                protein_sequence = gene_row.iloc[0]['protein_sequence']
                print(f"Found sequence for {gene_name} (Ensembl): {protein_sequence[:50]}...")
            else:
                protein_sequence = "None"
                print(f"No sequence found for gene '{gene_name}' in file '{file_path}' (Ensembl)")
            return protein_sequence
            

    except Exception as e:
        print(f"Error extracting sequence for {gene_name}: {e}")
        return None


def process_single_row(row, organism_file_map, mode="phytozome"):
    """
    Process a single row to extract gene sequence.
    
    :param row: DataFrame row
    :param organism_file_map: Mapping of organisms to file paths
    :param mode: 'phytozome' or 'ensembl'
    :return: Dictionary with extraction results
    """
    try:
        gene_name = row['gene_name']
        
        if row['organism'] in organism_file_map:
            file_path = organism_file_map[row['organism']]
            sequence = extract_gene_protein_sequence(file_path, gene_name, mode=mode)
            
            return {
                'gene_name': row['gene_name'],
                'organism': row['organism'],
                'protein_sequence': sequence
            }
        return None
    except Exception as e:
        print(f"Error processing row: {e}")
        return None


def extract_organism_file_mapping(csv_path, mode="phytozome"):
    """
    Create a dictionary mapping unique organism names to their corresponding file paths.
    
    :param csv_path: Path to the input CSV file
    :param mode: 'phytozome' or 'ensembl'
    :return: Dictionary of unique organism names to file paths
    """
    print("Starting organism file mapping...")
    start_time = time.time()
    
    if mode == "phytozome":
        base_path = "/groups/itay_mayrose/alongonda/full_genomes/phytozome/fasta_files/Phytozome"
    elif mode == "ensembl":
        base_path = "/groups/itay_mayrose/alongonda/full_genomes/ensembl/organisms"
    elif mode == "plaza":
        base_path = "/groups/itay_mayrose/alongonda/full_genomes/plaza/organisms"
    
    # Read CSV and get unique organisms
    df = pd.read_csv(csv_path, header=None, names=['gene_name', 'organism'])
    unique_organisms = df['organism'].unique()
    print(f"Total unique organisms: {len(unique_organisms)}")
    
    # Concurrent file mapping
    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        future_to_organism = {
            executor.submit(find_file_for_organism, base_path, org, mode): org 
            for org in unique_organisms
        }
        
        organism_file_map = {}
        for future in as_completed(future_to_organism):
            organism = future_to_organism[future]
            try:
                file_path = future.result()
                if file_path:
                    organism_file_map[organism] = file_path
            except Exception as exc:
                print(f"{organism} generated an exception: {exc}")
                break
    
    print(f"Organisms with files mapped: {len(organism_file_map)}")
    print(f"Mapping completed in {time.time() - start_time:.2f} seconds")
    
    return organism_file_map


def process_csv_and_extract_sequences(input_csv_path, output_fasta_path, mode="phytozome"):
    """
    Process CSV file to extract gene protein sequences using concurrent processing.
    
    :param input_csv_path: Path to input CSV file
    :param output_csv_path: Path to output CSV file for results
    :param mode: 'phytozome' or 'ensembl'
    """
    print(f"Starting sequence extraction from {input_csv_path} in mode {mode}")
    start_time = time.time()
    
    # Read input CSV
    df = pd.read_csv(input_csv_path, header=None, names=['gene_name', 'organism'])
    print(f"Total rows in input: {len(df)}")
    
    # Get organism to file mapping
    organism_file_map = extract_organism_file_mapping(input_csv_path, mode=mode)
    
    # Concurrent sequence extraction
    results = []
    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        futures = {
            executor.submit(process_single_row, row, organism_file_map, mode=mode): idx 
            for idx, row in df.iterrows()
        }
        
        for future in as_completed(futures):
            result = future.result()
            if result and result['protein_sequence']:
                results.append(result)
                
    with open(output_fasta_path, 'w') as fasta_file:
        for result in results:
            fasta_file.write(f'>{result["gene_name"]} | {result["organism"]}\n{result["protein_sequence"]}\n')
    
    print("\n--- Extraction Summary ---")
    print(f"Total rows processed: {len(df)}")
    print(f"Rows with sequences: {len(results)}")
    print(f"Sequences saved to: {output_fasta_path}")
    print(f"Total extraction time: {time.time() - start_time:.2f} seconds")


# Example usage
if __name__ == "__main__":
    input_csv_path = "/groups/itay_mayrose/alongonda/datasets/plantcyc/pmn_mgc_potential/phytozome/phytozome_genes.csv"  # Replace with your input CSV path
    output_fasta_path = "/groups/itay_mayrose/alongonda/datasets/plantcyc/pmn_mgc_potential/phytozome/phytozome_genes_with_sequence.fasta"  # Replace with desired output path
    mode = "phytozome"  # 'phytozome' or 'ensembl' or 'plaza'
    process_csv_and_extract_sequences(input_csv_path, output_fasta_path, mode=mode)
    print("Sequence extraction completed.")
