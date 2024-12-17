from threading import Lock
import os
import pandas as pd
import hashlib
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import ahocorasick
from tqdm import tqdm  # For progress bar

# Aho-Corasick Implementation to match multiple patterns in text
def build_aho_corasick(pathway_dict):
    automaton = ahocorasick.Automaton()
    for pathway, gene_lists in pathway_dict.items():
        for occurrence_idx, gene_list in enumerate(gene_lists):
            # Ensure the unique occurrence index is consistent
            unique_occurrence_id = f"{pathway}__{occurrence_idx}"
            for gene in gene_list:
                automaton.add_word(gene.lower(), (unique_occurrence_id, gene.lower()))
    automaton.make_automaton()
    return automaton

class UniqueMatchTracker:
    def __init__(self):
        self._seen_matches = []
        self._lock = Lock()

    def is_unique_match(self, pathway_cluster):
        new_cluster_genes = set(pathway_cluster['pathway_cluster_genes'].split(','))

        with self._lock:
            # Check if the new match is a subset of any existing match
            for match in self._seen_matches:
                existing_genes = set(match.split(','))
                if new_cluster_genes.issubset(existing_genes):
                    return False

            # If unique, add the new match
            self._seen_matches.append(pathway_cluster['pathway_cluster_genes'])
            return True

def generate_group_hash(group):
    sorted_genes = sorted(group['pathway_cluster_genes'].split(','))
    hash_string = f"{'|'.join(sorted_genes)}"
    return hashlib.md5(hash_string.encode()).hexdigest()

# Process a sliding window of data, matching gene sequences with Aho-Corasick
def process_window_with_aho_corasick(window_data, automaton, unique_tracker, output_file, file_lock, file_path, min_genes):
    gene_ids = [gene.lower() for gene in window_data['id']]
    matches_found = 0
    matched_pathways = defaultdict(lambda: defaultdict(set))  # Nested dict to track matches per occurrence
    gene_positions = defaultdict(dict)  # Store start and end positions for matched genes

    # Use Aho-Corasick to find all matching genes in the window
    for end_pos, (unique_occurrence_id, gene) in automaton.iter(" ".join(gene_ids)):
        pathway, occurrence_idx = unique_occurrence_id.split('__')
        matched_pathways[pathway][occurrence_idx].add(gene)

    # Process matches for each pathway occurrence
    for pathway, occurrences in matched_pathways.items():
        for occurrence_idx, genes in occurrences.items():
            # Find window genes that match the pathway (in their original case from window_data['id'])
            if len(genes) >= min_genes:
                matched_window_genes = []
                for window_gene, start, end in zip(window_data['id'], window_data['start'], window_data['end']):
                    if any(gene in window_gene.lower() for gene in genes):
                        matched_window_genes.append(window_gene)
                        gene_positions[window_gene] = {'start': start, 'end': end}

                unique_genes = set(genes)
                matched_positions = [
                    f"{gene}({gene_positions[gene]['start']}-{gene_positions[gene]['end']})"
                    for gene in matched_window_genes if gene in gene_positions
                ]
                group = {
                    'pathway': f"{pathway} (Occurrence {occurrence_idx})",  # Include occurrence info
                    'genes': ','.join(window_data['id']),
                    'pathway_cluster_genes': ','.join(sorted(unique_genes)),
                    'window_cluster_genes': ','.join(sorted(matched_window_genes)),  # Only matched genes from the window
                    'gene_positions': ','.join(matched_positions),  # Include gene positions
                    'start': window_data['start'].min(),
                    'end': window_data['end'].max(),
                    'source_file': file_path
                }
                pathway_cluster = {'pathway_cluster_genes': group['pathway_cluster_genes']}
                if unique_tracker.is_unique_match(pathway_cluster):
                    matches_found += 1
                    with file_lock:
                        mode = 'a' if os.path.exists(output_file) else 'w'
                        header = not os.path.exists(output_file)
                        pd.DataFrame([group]).to_csv(output_file, mode=mode, header=header, index=False)
                    print(f"Match Found: Pathway={pathway} (Occurrence {occurrence_idx}), "
                          f"Genes={group['genes']}, Matching Genes={group['pathway_cluster_genes']}, "
                          f"File={group['source_file']}")
    return matches_found

# Optimized create_pathways_dict
def create_pathways_dict(pathways_file):
    pathways_dict = defaultdict(list)  # Use a list to track each occurrence separately
    chunk_size = 10000  # Larger chunk size to optimize reading large files
    for chunk in pd.read_csv(pathways_file, chunksize=chunk_size, low_memory=False):
        # Filter relevant columns (columns with "GENE-ID." prefix and 'UNIQUE-ID')
        gene_id_cols = [col for col in chunk.columns if col.startswith('GENE-ID.')]

        # Prepare for batch processing: create a dict to map unique ids to gene ids
        for _, row in chunk.iterrows():
            unique_id = row['UNIQUE-ID']
            gene_ids = [str(val).lower() for col in gene_id_cols 
                        if pd.notna(val := row[col]) and val]
            if gene_ids:  # Only add to dict if there are valid gene IDs
                pathways_dict[unique_id].append(gene_ids)
    return dict(pathways_dict)

# Process each file in a directory using multithreading
def process_file(file_path, pathway_tries, output_file, unique_tracker, file_lock, window_size, min_genes=3):
    print(f"Processing file: {file_path} with window size: {window_size}")
    total_matches = 0
    df = pd.read_csv(file_path, usecols=['id', 'start', 'end'])
    df = df.sort_values(by=['start', 'end'])
    num_windows = len(df) - window_size + 1
    with tqdm(total=num_windows, desc=f"File: {os.path.basename(file_path)}", unit="window") as pbar:
        for i in range(num_windows):
            window = df.iloc[i:i + window_size]
            window_matches = process_window_with_aho_corasick(window, pathway_tries, unique_tracker, output_file, file_lock, file_path, min_genes)
            total_matches += window_matches
            pbar.update(1)
    print(f"Completed file: {file_path}, Matches Found: {total_matches}")
    return total_matches

# Process a whole genome directory
def process_genome_dir(genome_dir, pathway_tries, output_file, max_workers, window_size):
    file_paths = [os.path.join(genome_dir, f) for f in os.listdir(genome_dir) 
                  if f.endswith('.csv') and os.path.isfile(os.path.join(genome_dir, f))]
    unique_tracker = UniqueMatchTracker()
    file_lock = Lock()
    total_matches = 0
    with tqdm(total=len(file_paths), desc=f"Directory: {os.path.basename(genome_dir)}", unit="file") as pbar:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(
                    process_file, 
                    file_path, 
                    pathway_tries, 
                    output_file, 
                    unique_tracker, 
                    file_lock,
                    window_size
                )
                for file_path in file_paths
            ]
            for future in as_completed(futures):
                total_matches += future.result()
                pbar.update(1)
    print(f"Directory {genome_dir} - Total Matches: {total_matches}")
    return total_matches

def main():
    genome_dirs = [
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/ensembl/processed_annotations_sorted",
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/plaza/processed_annotations_sorted",
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/phytozome/processed_annotations"
    ]
    pathways_file = "/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms/merged_pathways.csv"
    output_dir = "/groups/itay_mayrose_nosnap/alongonda/Plant_MGC/sliding_window_outputs"
    
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    max_workers = min(32, os.cpu_count())
    print(f"Using {max_workers} workers")
    print("Loading pathways...")

    # Load pathway dictionary and build Aho-Corasick automaton
    pathway_dict = create_pathways_dict(pathways_file)
    pathway_tries = build_aho_corasick(pathway_dict)
    print(f"Loaded {len(pathway_dict)} pathways with Aho-Corasick.")

    for window_size in range(18, 21):
        output_file = os.path.join(output_dir, f"potential_groups_with_window_size_{window_size}.csv")
        if os.path.exists(output_file):
            os.remove(output_file)

        total_matches = 0
        for genome_dir in genome_dirs:
            print(f"Processing directory: {genome_dir} with window size: {window_size}")
            total_matches += process_genome_dir(genome_dir, pathway_tries, output_file, max_workers, window_size)
        
        print(f"TOTAL MATCHES FOUND for window size {window_size}: {total_matches}")
        print(f"Results saved to: {output_file}")

if __name__ == "__main__":
    main()