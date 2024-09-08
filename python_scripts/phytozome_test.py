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
        print(unique_id)
        query = service.new_query("Protein")
        query.add_view("sequence.residues")
        query.add_constraint("primaryIdentifier", "=", unique_id, code="A")
        for row in query.rows():
            return f'{unique_id} (gene {gene_to_id[unique_id]}): {row["sequence.residues"]}\n'
    except Exception as e:
        print(f"Error processing unique_id {unique_id}: {e}")
        return None

def search_in_phytozome_and_write_to_file(unique_ids, gene_to_id, unique_ids_seq):
    service = Service("https://phytozome-next.jgi.doe.gov/phytomine/service")
    
    with open(unique_ids_seq, 'w') as unique_ids_seq_file:
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures = {executor.submit(process_unique_id, service, unique_id, gene_to_id): unique_id for unique_id in unique_ids}
            for future in as_completed(futures):
                result = future.result()
                if result:
                    unique_ids_seq_file.write(result)

def main():
    unique_ids_file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/gene_by_db_seq/unique_dbs_for_unique_id.txt'
    unique_ids_seq = '/groups/itay_mayrose/alongonda/desktop/plantcyc/gene_by_db_seq/PHYTOZOME/unique_ids_seq.txt'
    
    # Extract UNIQUE-IDs with ENSEMBL UNIQUE-DBLINK
    unique_ids, gene_to_id = extract_unique_ids(unique_ids_file_path)
    
    # Search for UNIQUE-IDs in Phytozome and write to file
    search_in_phytozome_and_write_to_file(unique_ids, gene_to_id, unique_ids_seq)

if __name__ == '__main__':
    main()