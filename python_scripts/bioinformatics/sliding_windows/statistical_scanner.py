from threading import Lock
import os
import numpy as np
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

def clean_and_process_genes_and_pathways(genes_and_pathways):
    processed = {}
    for gene, values in genes_and_pathways.items():
        # Remove unwanted characters and split into elements
        stripped_values = [''.join(c for c in v if c not in "[]'\"").split(",") for v in values]
        processed[gene] = stripped_values
    return processed

def find_first_common_element(genes_and_pathways, min_genes):
    from itertools import product, combinations
    # Clean and process the input
    processed = clean_and_process_genes_and_pathways(genes_and_pathways)
    
    # Flatten any nested lists within each gene's pathways
    sets = []
    for paths in processed.values():
        flattened = set(item.strip() for sublist in paths for item in sublist)
        sets.append(flattened)
        
    # Generate all combinations of sets of size min_genes
    set_combinations = combinations(sets, min_genes)
    
    for selected_sets in set_combinations:
        # Generate all combinations with one element per selected set
        all_combinations = product(*selected_sets)
        
        for combination in all_combinations:
            # If all elements are the same, return the first match
            if len(set(combination)) == 1:
                return combination[0]
    
    return None  # No common element found

def process_df(filtered_df, window_size, min_genes, file_path, file_lock, output_file):
    total_matches = 0
    # Group by chromosome and process each group
    chromosomes = filtered_df['chromosome'].unique()
    # Process each chromosome group
    for chromosome in chromosomes:
        chromosome_data = filtered_df[filtered_df['chromosome'] == chromosome]
        num_genes = len(chromosome_data)
        i = 0
        while i < num_genes:
            genes_and_pathways = dict()
            window = [chromosome_data.iloc[i]]
            genes_and_pathways[chromosome_data.iloc[i]['metabolic_gene']] = [chromosome_data.iloc[i]['pathway']]
            start_index = chromosome_data.iloc[i]['index']
            for j in range(i+1, num_genes):
                end_index = chromosome_data.iloc[j]['index']
                if (end_index - start_index <= window_size):
                    window.append(chromosome_data.iloc[j])
                    if (chromosome_data.iloc[j]['metabolic_gene'] not in genes_and_pathways.keys()):
                        genes_and_pathways[chromosome_data.iloc[j]['metabolic_gene']] = [chromosome_data.iloc[j]['pathway']]
                    else:
                        genes_and_pathways[chromosome_data.iloc[j]['metabolic_gene']].append(chromosome_data.iloc[j]['pathway'])
                else:
                    break
            if len(genes_and_pathways.keys()) >= min_genes:
                window_df = pd.DataFrame(window)
                pathway = find_first_common_element(genes_and_pathways, min_genes)
                if pathway:
                    group = {
                        'pathway': pathway,  # Include occurrence info
                        'genes': ','.join(window_df['id']),
                        'metabolic_genes': ','.join(list(genes_and_pathways.keys())),
                        'start': window_df['start'].min(),
                        'end': window_df['end'].max(),
                        'source_file': file_path
                    }
                    with file_lock:
                        mode = 'a' if os.path.exists(output_file) else 'w'
                        header = not os.path.exists(output_file)
                        pd.DataFrame([group]).to_csv(output_file, mode=mode, header=header, index=False)
                        total_matches += 1
            i += 1
    return total_matches

# Process each file in a directory using multithreading
def process_file(file_path, output_file, file_lock, window_size, min_genes):
    print(f"Processing file: {file_path} with window size: {window_size}")
    df = pd.read_csv(file_path)
    df["index"] = df.index
    # Filter out non-metabolic genes
    filtered_df = df[df["pathway"].notna()]
    total_matches = process_df(filtered_df, window_size, min_genes, file_path, file_lock, output_file)
    
    random_advantage = 10
    # Shuffle the values in columns 'pathway' and 'metabolic_gene' identically
    shuffled_df = df.copy()
    total_random_matches = 0
    for i in range(random_advantage):
        permutation = np.random.permutation(len(shuffled_df))
        shuffled_df['pathway'] = shuffled_df['pathway'].values[permutation]
        shuffled_df['metabolic_gene'] = shuffled_df['metabolic_gene'].values[permutation]
        filtered_shuffled_df = shuffled_df[shuffled_df["pathway"].notna()]
        
        total_random_matches += process_df(filtered_shuffled_df, window_size, min_genes, file_path, file_lock, output_file) / random_advantage
    print(f"Completed file: {file_path}, Matches Found: {total_matches}")
    return total_matches, total_random_matches

# Process a whole genome directory
def process_genome_dir(genome_dir, output_file, max_workers, window_size, min_genes):
    file_paths = [os.path.join(genome_dir, f) for f in os.listdir(genome_dir) 
                  if f.endswith('.csv') and os.path.isfile(os.path.join(genome_dir, f))]
    file_lock = Lock()
    total_matches = 0
    total_random_matches = 0
    with tqdm(total=len(file_paths), desc=f"Directory: {os.path.basename(genome_dir)}", unit="file") as pbar:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(
                    process_file, 
                    file_path, 
                    output_file, 
                    file_lock,
                    window_size,
                    min_genes
                )
                for file_path in file_paths
            ]
            for future in as_completed(futures):
                future_matches, future_random_matches = future.result()
                total_matches += future_matches
                total_random_matches += future_random_matches
                pbar.update(1)
    print(f"Directory {genome_dir} - Total Matches: {total_matches}")
    print(f"Directory {genome_dir} - Total Random Matches: {total_random_matches}")
    return total_matches, total_random_matches

def create_output_subdir(output_dir, min_genes):
    """Create a subdirectory for each min_genes value."""
    subdir = os.path.join(output_dir, f"min_genes_{min_genes}")
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    return subdir

def enrichment_analysis(total_random_matches, total_matches, window_size, min_genes, enrichment_output_file):
    # Calculate enrichment statistics
    if total_random_matches > 0:
        enrichment_ratio = total_matches / total_random_matches
    else:
        enrichment_ratio = float('inf')  # Avoid division by zero

    enrichment_stats = {
        'window_size': window_size,
        'min_genes': min_genes,
        'total_matches': total_matches,
        'total_random_matches': total_random_matches,
        'enrichment_ratio': enrichment_ratio
    }

    # Save enrichment statistics to a CSV file
    enrichment_df = pd.DataFrame([enrichment_stats])
    mode = 'a' if os.path.exists(enrichment_output_file) else 'w'
    header = not os.path.exists(enrichment_output_file)
    enrichment_df.to_csv(enrichment_output_file, mode=mode, header=header, index=False)

def main():
    # genome_dirs = [
    #     "/groups/itay_mayrose/alongonda/full_genomes/ensembl/processed_annotations_test_no_chloroplast_with_sequences",
    #     "/groups/itay_mayrose/alongonda/full_genomes/plaza/processed_annotations_with_chromosomes_no_chloroplast_with_sequences",
    #     "/groups/itay_mayrose/alongonda/full_genomes/phytozome/processed_annotations_with_chromosomes_no_chloroplast_with_sequences"
    # ]
    genome_dirs = ["/groups/itay_mayrose/alongonda/full_genomes/mgc_enriched_files"]
    pathways_file = "/groups/itay_mayrose/alongonda/plantcyc/all_organisms/merged_pathways.csv"
    output_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/10_most_enriched_genomes_statistical_tests"
    
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    max_workers = min(32, os.cpu_count())
    print(f"Using {max_workers} workers")

    for window_size in range(5, 21):
        # Dynamically calculate the maximum value of min_genes
        max_min_genes = (window_size // 2) + 1
        
        # for min_genes in range(3, max_min_genes + 1):  # Adjust min_genes based on the constraint
        for min_genes in range(3, max_min_genes):
            # Create subdirectory for the current min_genes
            min_genes_subdir = create_output_subdir(output_dir, min_genes)
            output_file = os.path.join(min_genes_subdir, f"potential_groups_w{window_size}g{min_genes}.csv")
            enrichment_output_file = os.path.join(min_genes_subdir, f"enrichment_w{window_size}g{min_genes}.csv")
            
            if os.path.exists(output_file):
                os.remove(output_file)

            total_matches = 0
            total_random_matches = 0
            for genome_dir in genome_dirs:
                print(f"Processing genome directory: {genome_dir} with window size: {window_size} and min_genes: {min_genes}")
                matches, random_matches = process_genome_dir(genome_dir, output_file, max_workers, window_size, min_genes)
                total_matches += matches
                total_random_matches += random_matches
            
            print(f"TOTAL MATCHES FOUND for window size {window_size} and min_genes {min_genes}: {total_matches}")
            print(f"TOTAL RANDOM MATCHES FOUND for window size {window_size} and min_genes {min_genes}: {total_random_matches}")
            
            enrichment_analysis(total_random_matches, total_matches, window_size, min_genes, enrichment_output_file)
            
            print(f"Results saved to: {output_file}")
            print(f"Enrichment statistics saved to: {enrichment_output_file}")
            
if __name__ == "__main__":
    main()