import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

def extract_unique_ids(file_path):
    unique_ids = set()
    gene_to_id = {}
    try:
        with open(file_path, 'r', errors='ignore') as file:
            lines = file.readlines()
        for line in lines:
            if 'ENSEMBL' in line:
                match = re.search(r'UNIQUE-ID: (\S+)', line)
                gene = re.search(r'GENE: (\S+)', line)
                if match:
                    value = match.group(1).replace('-MONOMER','').replace('MONOMER','').replace(',','')
                    unique_ids.add(value)
                    gene_to_id[value] = gene.group(1).replace(',', '')
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
    return unique_ids, gene_to_id

def process_fasta_file(file_path, unique_ids):
    ids_and_sequences = {}
    found_ids = set()
    try:
        with open(file_path, 'r', errors='ignore') as fasta_file:
            print(f"Processing file: {file_path}")
            current_id = None
            for line in fasta_file:
                if line.startswith('>transcript:'):
                    current_id = line.split(':')[1].strip()
                elif current_id and current_id in unique_ids:
                    ids_and_sequences[current_id] = line.strip() if current_id not in ids_and_sequences else f"{ids_and_sequences[current_id]}{line.strip()}"
                    found_ids.add(current_id)
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
    return ids_and_sequences, found_ids

def search_in_fasta_files(unique_ids, base_dir, gene_to_id):
    ids_and_sequences_general = {}
    all_found_ids = set()
    file_paths = []

    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith('.canon.fa'):
                file_paths.append(os.path.join(root, file))

    with ThreadPoolExecutor(max_workers=16) as executor:
        futures = {executor.submit(process_fasta_file, file_path, unique_ids): file_path for file_path in file_paths}
        for future in as_completed(futures):
            ids_and_sequences, found_ids = future.result()
            ids_and_sequences_general.update(ids_and_sequences)
            all_found_ids.update(found_ids)

    unique_ids -= all_found_ids  # Remove all found unique_ids at once
    
    genes_and_sequences = {}
    for id, sequence in ids_and_sequences_general.items():
        gene_id = gene_to_id[id]
        if (gene_id in genes_and_sequences and len(genes_and_sequences[gene_id]) < len(sequence)) or gene_id not in genes_and_sequences:
            genes_and_sequences[gene_id] = sequence
    return genes_and_sequences, unique_ids

def main():
    unique_ids_file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/gene_by_db_seq/unique_dbs_for_unique_id.txt'
    fasta_base_dir = '/groups/itay_mayrose/galtoledano/gene_families/data/snakemake_download/per_species'
    unique_ids_seq = '/groups/itay_mayrose/alongonda/desktop/plantcyc/gene_by_db_seq/ENSEMBL/unique_ids_seq.txt'
    
    # Extract UNIQUE-IDs with ENSEMBL UNIQUE-DBLINK
    unique_ids, gene_to_id = extract_unique_ids(unique_ids_file_path)
    
    # Search for UNIQUE-IDs in .canon.fa files
    genes_and_sequences, unfound_ids = search_in_fasta_files(unique_ids, fasta_base_dir, gene_to_id)
    # Write the results to the output file
    with open(unique_ids_seq, 'w') as unique_ids_seq_file:
        for gene, sequence in genes_and_sequences.items():
            unique_ids_seq_file.write(f'{gene}: {sequence}\n')
            
        for unfound_id in unfound_ids:
            unique_ids_seq_file.write(f'{gene_to_id[unfound_id]}: N/A\n')

if __name__ == '__main__':
    main()