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
    output_file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms/dblinks_for_uniqueid.txt'
    
    results = []

    # Walk through all directories and files under base_dir
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file == 'proteins.dat':
                file_path = os.path.join(root, file)
                sub_db = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(file_path))))
                print(f"Processing file: {file_path} for sub_db: {sub_db}")
                matches = extract_dblinks(file_path)
                if matches:
                    for unique_id, dblinks in matches:
                        results.append((sub_db, unique_id, dblinks))

    # Write the results to the output file
    with open(output_file_path, 'w') as output_file:
        for sub_db, unique_id, dblinks in results:
            output_file.write(f'sub_db: {sub_db}, UNIQUE-ID: {unique_id}, DBLINKS: {dblinks}\n')

if __name__ == '__main__':
    main()