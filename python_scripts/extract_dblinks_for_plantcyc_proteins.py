import os
import re

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
            print(f"No matches found in file: {file_path}")

        return matches
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return []

def main():
    base_dir = '/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms'
    output_file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms/dblinks_for_uniqueids.txt'
    unique_dbs_file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms/unique_dbs_count.txt'
    unique_dbs_for_unique_id_file_path = '/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms/unique_dbs_for_unique_id.txt'
    
    results = []
    unique_dbs = {}

    # Walk through all directories and files under base_dir
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file == 'proteins.dat':
                file_path = os.path.join(root, file)
                organism = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(file_path))))
                print(f"Processing file: {file_path} for organism: {organism}")
                matches = extract_dblinks(file_path)
                if matches:
                    for unique_id, gene, dblinks in matches:
                        results.append((organism, unique_id, gene, dblinks))
                        # Extract the database name and add to the set of unique DBs
                        db_names = re.findall(r'\((\w+)', dblinks)
                        for db_name in db_names:
                            if db_name in unique_dbs:
                                unique_dbs[db_name] += 1
                            else:
                                unique_dbs[db_name] = 1

    # Sort unique_dbs by count in descending order
    sorted_unique_dbs = sorted(unique_dbs.items(), key=lambda item: item[1], reverse=True)
    
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
                        unique_dbs[db_name] -= 1
                uniques.append((organism, unique_id, gene, unique_db))
            else:
                new_results.append((organism, unique_id, gene, dblinks))
        results_copy = new_results

        # Re-sort unique_dbs after processing each unique_db
        sorted_unique_dbs = sorted(unique_dbs.items(), key=lambda item: item[1], reverse=True)
        index += 1

    # Write the results to the output file
    with open(output_file_path, 'w') as output_file:
        for organism, unique_id, gene, dblinks in results:
            output_file.write(f'Organism: {organism}, UNIQUE-ID: {unique_id}, GENE: {gene}, DBLINKS: {dblinks}\n')

    # Write the unique database names to the unique_dbs file
    with open(unique_dbs_file_path, 'w') as unique_dbs_file:
        for db_name, count in sorted(unique_dbs.items(), key=lambda item: item[1], reverse=True):
            unique_dbs_file.write(f'{db_name}: {count}\n')
            
    # Write the results to the output file
    with open(unique_dbs_for_unique_id_file_path, 'w') as unique_dbs_for_unique_id:
        for organism, unique_id, gene, unique_db in uniques:
            unique_dbs_for_unique_id.write(f'Organism: {organism}, UNIQUE-ID: {unique_id}, GENE: {gene}, UNIQUE-DBLINK: {unique_db}\n')

if __name__ == '__main__':
    main()