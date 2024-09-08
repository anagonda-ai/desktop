import os
import re

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

def search_in_fasta_files(unique_ids, base_dir):
    results = {}
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith('.canon.fa'):
                file_path = os.path.join(root, file)
                try:
                    with open(file_path, 'r', errors='ignore') as fasta_file:
                        print(f"Processing file: {file_path}")
                        found_ids = set()
                        current_id = None
                        for line in fasta_file:
                            if line.startswith('>transcript:'):
                                current_id = line.split(':')[1].strip()
                            elif current_id and current_id in unique_ids:
                                results[current_id] = line.strip() if current_id not in results else f"{results[current_id]}{line.strip()}"
                                found_ids.add(current_id)
                        unique_ids -= found_ids  # Remove all found unique_ids at once
                except Exception as e:
                    print(f"Error processing file {file_path}: {e}")
    return results, unique_ids

def main():
    unique_ids_file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/gene_by_db_seq/unique_dbs_for_unique_id.txt'
    fasta_base_dir = '/groups/itay_mayrose/galtoledano/gene_families/data/snakemake_download/per_species'
    unique_ids_seq = '/groups/itay_mayrose/alongonda/desktop/plantcyc/gene_by_db_seq/ENSEMBL/unique_ids_seq.txt'
    
    # Extract UNIQUE-IDs with ENSEMBL UNIQUE-DBLINK
    unique_ids, gene_to_id = extract_unique_ids(unique_ids_file_path)
    
    # Search for UNIQUE-IDs in .canon.fa files
    results, unfounds = search_in_fasta_files(unique_ids, fasta_base_dir)
    # Write the results to the output file
    with open(unique_ids_seq, 'w') as unique_ids_seq_file:
        for unique_id, value in results.items():
            unique_ids_seq_file.write(f'{unique_id} (gene {gene_to_id[unique_id]}): {value}\n')
            
        for unfound in unfounds:
            unique_ids_seq_file.write(f'{unfound} (gene {gene_to_id[unique_id]}): N\A\n')

if __name__ == '__main__':
    main()