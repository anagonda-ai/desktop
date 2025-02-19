from threading import Lock
import os
import pandas as pd
import hashlib
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm  # For progress bar

def generate_group_hash(group):
    sorted_genes = sorted(group['pathway_cluster_genes'].split(','))
    hash_string = f"{'|'.join(sorted_genes)}"
    return hashlib.md5(hash_string.encode()).hexdigest()


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
    gene_sets = {}
    for gene, paths in processed.items():
        flattened = set(item.strip() for sublist in paths for item in sublist)
        gene_sets[gene] = flattened
    
    # Generate all combinations of genes of size min_genes
    gene_combinations = combinations(gene_sets.keys(), min_genes)
    
    for selected_genes in gene_combinations:
        selected_sets = [gene_sets[gene] for gene in selected_genes]
        
        # Generate all combinations with one element per selected set
        all_combinations = product(*selected_sets)
        
        for combination in all_combinations:
            # If all elements are the same, return the first match and corresponding genes
            if len(set(combination)) == 1:
                return combination[0], list(selected_genes)
    
    return None, []  # No common element found

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
def process_file(file_path, output_file, file_lock, window_size, min_genes):
    print(f"Processing file: {file_path} with window size: {window_size}")
    total_matches = 0
    df = pd.read_csv(file_path)
    df["index"] = df.index
    # Filter out non-metabolic genes
    filtered_df = df[df["pathway"].notna()]
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
                pathway, metabolic_genes = find_first_common_element(genes_and_pathways, min_genes)
                if pathway:
                    group = {
                        'pathway': pathway,  # Include occurrence info
                        'genes': ','.join(window_df['id']),
                        'metabolic_genes': ','.join(metabolic_genes),
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
    print(f"Completed file: {file_path}, Matches Found: {total_matches}")
    return total_matches

# Process a whole genome directory
def process_genome_dir(genome_dir, output_file, max_workers, window_size, min_genes):
    file_paths = [os.path.join(genome_dir, f) for f in os.listdir(genome_dir) 
                  if f.endswith('.csv') and os.path.isfile(os.path.join(genome_dir, f))]
    file_lock = Lock()
    total_matches = 0
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
                total_matches += future.result()
                pbar.update(1)
    print(f"Directory {genome_dir} - Total Matches: {total_matches}")
    return total_matches

def create_output_subdir(output_dir, min_genes):
    """Create a subdirectory for each min_genes value."""
    subdir = os.path.join(output_dir, f"min_genes_{min_genes}")
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    return subdir

def main():
    genome_dirs = [
        "/groups/itay_mayrose/alongonda/full_genomes/ensembl/processed_annotations_test_no_chloroplast_with_sequences/metabolic_genes_with_letters",
        "/groups/itay_mayrose/alongonda/full_genomes/plaza/processed_annotations_with_chromosomes_no_chloroplast_with_sequences/metabolic_genes_with_letters",
        "/groups/itay_mayrose/alongonda/full_genomes/phytozome/processed_annotations_with_chromosomes_no_chloroplast_with_sequences/metabolic_genes_with_letters"
    ]
    pathways_file = "/groups/itay_mayrose/alongonda/datasets/plantcyc/all_organisms/merged_pathways.csv"
    output_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/min_genes_only_metabolic_genes_input_test"
    
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    max_workers = min(32, os.cpu_count())
    print(f"Using {max_workers} workers")
    print("Loading pathways...")

    for window_size in range(5, 21):
        # Dynamically calculate the maximum value of min_genes
        max_min_genes = (window_size // 2) + 1
        
        # for min_genes in range(3, max_min_genes + 1):  # Adjust min_genes based on the constraint
        for min_genes in [3]:
            # Create subdirectory for the current min_genes
            min_genes_subdir = create_output_subdir(output_dir, min_genes)
            output_file = os.path.join(min_genes_subdir, f"potential_groups_w{window_size}.csv")
            
            if os.path.exists(output_file):
                os.remove(output_file)

            total_matches = 0
            for genome_dir in genome_dirs:
                print(f"Processing genome directory: {genome_dir} with window size: {window_size} and min_genes: {min_genes}")
                process_genome_dir(genome_dir, output_file, max_workers, window_size, min_genes)
            
            print(f"TOTAL MATCHES FOUND for window size {window_size} and min_genes {min_genes}: {total_matches}")
            print(f"Results saved to: {output_file}")
            
if __name__ == "__main__":
    main()