
class DblinkExtractorProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: dblink_extractor.py."""
    
    def __init__(self, **kwargs):
        """Initialize processor."""
        super().__init__(**kwargs)
        self.settings = get_settings()
    
    def validate_input(self, data):
        """Validate input data."""
        pass  # Implement validation
    
    def process(self, data, **kwargs):
        """Process data."""
        # Original script logic here
        pass


from ..config.settings import get_settings
from ..core.base import BioinformaticsProcessor
from ..core.exceptions import ProcessingError, ValidationError
from ..utils.file_operations import file_manager
from loguru import logger
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

def extract_dblinks(file_path):
    try:
        with open(file_path, 'r', errors='ignore') as file:
            lines = file.readlines()

        matches = []
        unique_id = None
        dblinks_list = []
        gene = None
        for line in lines:
            if line.startswith('UNIQUE-ID - '):
                if unique_id and dblinks_list:
                    # Append all collected DBLINKS for the previous UNIQUE-ID
                    matches.append((unique_id, gene, ' '.join(dblinks_list)))
                unique_id = line.split(' - ')[1].strip()
                dblinks_list = []  # Reset the DBLINKS list for the new UNIQUE-ID
            elif line.startswith('DBLINKS - ') and unique_id:
                dblinks = line.split(' - ')[1].strip()
                dblinks_list.append(dblinks)
            elif line.startswith('GENE - ') and unique_id:
                gene = line.split(' - ')[1].strip()

        # Append the last collected DBLINKS if any
        if unique_id and dblinks_list:
            matches.append((unique_id, gene, ' '.join(dblinks_list)))

        if not matches:
            logger.info(f"No matches found in file: {file_path}")

        return matches
    except Exception as e:
        logger.info(f"Error processing file {file_path}: {e}")
        return []

def process_file(file_path):
    organism = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(file_path))))
    logger.info(f"Processing file: {file_path} for organism: {organism}")
    matches = extract_dblinks(file_path)
    return organism, matches

def main():
    base_dir = '/groups/itay_mayrose/alongonda/datasets/plantcyc/all_organisms'
    output_file_path = '/groups/itay_mayrose/alongonda/datasets/plantcyc/gene_by_db_seq/dblinks_for_uniqueids.txt'
    unique_dbs_file_path = '/groups/itay_mayrose/alongonda/datasets/plantcyc/gene_by_db_seq/unique_dbs_count.txt'
    unique_dbs_for_unique_id_file_path = '/groups/itay_mayrose/alongonda/datasets/plantcyc/gene_by_db_seq/unique_dbs_for_unique_id.txt'
    
    results = []
    unique_dbs = {}

    # Collect all file paths to process
    file_paths = []
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file == 'proteins.dat':
                file_paths.append(os.path.join(root, file))

    # Use ThreadPoolExecutor to process files in parallel
    with ThreadPoolExecutor(max_workers=4) as executor:
        future_to_file = {executor.submit(process_file, file_path): file_path for file_path in file_paths}
        for future in as_completed(future_to_file):
            file_path = future_to_file[future]
            try:
                organism, matches = future.result()
                if matches:
                    for unique_id, gene, dblinks in matches:
                        results.append((organism, unique_id, gene, dblinks))
                        # Extract the database name and add to the set of unique DBs
                        db_names = re.findall(r'\((\w+)', dblinks)
                        for db_name in db_names:
                            if db_name in unique_dbs:
                                unique_dbs[db_name].add(gene) # Add the unique_id without the transcript suffix. Just the protein
                            else:
                                unique_dbs[db_name] = set()
                                unique_dbs[db_name].add(gene)
            except Exception as e:
                logger.info(f"Error processing file {file_path}: {e}")

    # Sort unique_dbs by count in descending order
    sorted_unique_dbs = sorted(unique_dbs.items(), key=lambda item: len(item[1]), reverse=True)
    
    # Create copy of results to iterate over and remove DBLINKS that are not unique
    results_copy = results.copy()
    uniques = []
    
    # Iterate through each unique_db in sorted order
    index = 0
    while index < len(sorted_unique_dbs):
        unique_db, _ = sorted_unique_dbs[index]
        new_results = []
        for organism, unique_id, gene, dblinks in results_copy:
            db_names = re.findall(r'\((\w+)', dblinks)
            if unique_db in db_names:
                for db_name in db_names:
                    if db_name != unique_db:
                        unique_dbs[db_name].discard(gene)
                uniques.append((organism, unique_id, gene, unique_db))
            else:
                new_results.append((organism, unique_id, gene, dblinks))
        results_copy = new_results

        # Re-sort unique_dbs after processing each unique_db
        sorted_unique_dbs = sorted(unique_dbs.items(), key=lambda item: len(item[1]), reverse=True)
        index += 1

    # Write the results to the output file
    with open(output_file_path, 'w') as output_file:
        for organism, unique_id, gene, dblinks in results:
            output_file.write(f'Organism: {organism}, UNIQUE-ID: {unique_id}, GENE: {gene}, DBLINKS: {dblinks}\n')

    # Write the unique database names to the unique_dbs file
    with open(unique_dbs_file_path, 'w') as unique_dbs_file:
        for db_name, id_set in sorted(unique_dbs.items(), key=lambda item: len(item[1]), reverse=True):
            unique_dbs_file.write(f'{db_name}: {len(id_set)}\n')
            
    # Write the results to the output file
    with open(unique_dbs_for_unique_id_file_path, 'w') as unique_dbs_for_unique_id:
        for organism, unique_id, gene, unique_db in uniques:
            unique_dbs_for_unique_id.write(f'Organism: {organism}, UNIQUE-ID: {unique_id}, GENE: {gene}, UNIQUE-DBLINK: {unique_db}\n')

if __name__ == '__main__':
    main()