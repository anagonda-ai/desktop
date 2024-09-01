import os
import re

def extract_dblinks(file_path):
    try:
        with open(file_path, 'r', errors='ignore') as file:
            lines = file.readlines()

        matches = []
        unique_id = None
        for line in lines:
            if line.startswith('UNIQUE-ID - '):
                unique_id = line.split(' - ')[1].strip()
            elif line.startswith('DBLINKS - ') and unique_id:
                dblinks = line.split(' - ')[1].strip()
                matches.append((unique_id, dblinks))
                unique_id = None  # Reset unique_id after finding the first DBLINKS

        if not matches:
            print(f"No matches found in file: {file_path}")

        return matches
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return []

def main():
    base_dir = '/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms'
    output_file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms/dblinks_for_uniqueids.txt'
    unique_dbs_file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms/unique_dbs.txt'
    
    results = []
    unique_dbs = set()

    # Walk through all directories and files under base_dir
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file == 'proteins.dat':
                file_path = os.path.join(root, file)
                organism = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(file_path))))
                print(f"Processing file: {file_path} for organism: {organism}")
                matches = extract_dblinks(file_path)
                if matches:
                    for unique_id, dblinks in matches:
                        results.append((organism, unique_id, dblinks))
                        # Extract the database name and add to the set of unique DBs
                        match = re.match(r'\((\w+)', dblinks)
                        if match:
                            db_name = match.group(1)
                            unique_dbs.add(db_name)

    # Write the results to the output file
    with open(output_file_path, 'w') as output_file:
        for organism, unique_id, dblinks in results:
            output_file.write(f'Organism: {organism}, UNIQUE-ID: {unique_id}, DBLINKS: {dblinks}\n')

    # Write the unique database names to the unique_dbs file
    with open(unique_dbs_file_path, 'w') as unique_dbs_file:
        for db_name in sorted(unique_dbs):
            unique_dbs_file.write(f'{db_name}\n')

if __name__ == '__main__':
    main()