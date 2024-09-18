import os
import re
from intermine.webservice import Service
from concurrent.futures import ThreadPoolExecutor, as_completed

def extract_unique_ids(file_path):
    unique_ids = set()
    gene_to_id = {}
    try:
        with open(file_path, 'r', errors='ignore') as file:
            lines = file.readlines()
        for line in lines:
            if 'PHYTOZOME' in line:
                match = re.search(r'UNIQUE-ID: (\S+)', line)
                gene = re.search(r'GENE: (\S+)', line)
                if match:
                    value = match.group(1).replace('-MONOMER','').replace('MONOMER','').replace(',','')
                    unique_ids.add(value)
                    gene_to_id[value] = gene.group(1).replace(',', '')
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
    return unique_ids, gene_to_id

def process_unique_id(service, unique_id, gene_to_id):
    try:
        query = service.new_query("Protein")
        query.add_view("sequence.residues")
        query.add_constraint("primaryIdentifier", "=", unique_id, code="A")
        for row in query.rows():
            gene = gene_to_id[unique_id]
            sequence = row["sequence.residues"]
            return gene, sequence
    except Exception as e:
        print(f"Error processing unique_id {unique_id}: {e}")

def search_in_phytozome_and_write_to_file(unique_ids, gene_to_id, unique_ids_seq, no_duplicates):
    service = Service("https://phytozome-next.jgi.doe.gov/phytomine/service")
    
    genes_and_sequences = {}
    with open(unique_ids_seq, 'w') as unique_ids_seq_file:
        with ThreadPoolExecutor() as executor:
            futures = {executor.submit(process_unique_id, service, unique_id, gene_to_id): unique_id for unique_id in unique_ids}
            for future in as_completed(futures):
                result = future.result()
                print(f"Result: {result}")
                if result:
                    gene, sequence = result
                    if not no_duplicates:
                        unique_ids_seq_file.write(f">{gene}\n{sequence}\n")
                    else:
                        if gene and ((gene in genes_and_sequences and len(genes_and_sequences[gene]) < len(sequence)) or gene not in genes_and_sequences):
                                genes_and_sequences[gene] = sequence
                        
            if no_duplicates:
                for gene, sequence in genes_and_sequences.items():
                    unique_ids_seq_file.write(f">{gene}\n{sequence}\n")

def main():
    unique_ids_file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/gene_by_db_seq/unique_dbs_for_unique_id.txt'
    unique_ids_seq = '/groups/itay_mayrose/alongonda/desktop/plantcyc/gene_by_db_seq/PHYTOZOME/unique_ids_seq_test.txt'
    
    # Extract UNIQUE-IDs with ENSEMBL UNIQUE-DBLINK
    unique_ids, gene_to_id = extract_unique_ids(unique_ids_file_path)
    
    no_duplicates = False
    
    # Search for UNIQUE-IDs in Phytozome and write to file
    search_in_phytozome_and_write_to_file(unique_ids, gene_to_id, unique_ids_seq, no_duplicates)

if __name__ == '__main__':
    main()