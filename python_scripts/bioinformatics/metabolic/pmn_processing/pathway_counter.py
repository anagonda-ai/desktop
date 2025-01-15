import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

def count_unique_ids(file_path):
    unique_ids = defaultdict(int)
    try:
        with open(file_path, 'r', errors='ignore') as file:
            lines = file.readlines()
        for line in lines:
            if line.strip() and not (line.startswith('#') or line.startswith('UNIQUE-ID')):
                unique_id = line.split('\t')[0]
                unique_ids[unique_id] += 1
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
    return unique_ids

def main():
    base_dir = '/groups/itay_mayrose/alongonda/plantcyc/all_organisms'
    pathway_counts = defaultdict(int)
    file_paths = []

    # Collect all pathways.col file paths
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file == 'pathways.col':
                file_paths.append(os.path.join(root, file))

    # Use ThreadPoolExecutor to process files in parallel
    with ThreadPoolExecutor(max_workers=16) as executor:
        futures = {executor.submit(count_unique_ids, file_path): file_path for file_path in file_paths}
        for future in as_completed(futures):
            file_path = futures[future]
            try:
                unique_ids = future.result()
                for unique_id, count in unique_ids.items():
                    pathway_counts[unique_id] += count
                print(f"Processed {file_path}, found {len(unique_ids)} unique IDs")
            except Exception as e:
                print(f"Error processing file {file_path}: {e}")
    
    print(f"Total unique pathways found: {len(pathway_counts)}")

if __name__ == '__main__':
    main()