from threading import Lock
import string
import os
import pandas as pd
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm  # For progress bar

pathway_occurences = dict()

def increase_unique_id_occurence(unique_id):
    if unique_id not in pathway_occurences:
        pathway_occurences[unique_id] = 0
    pathway_occurences[unique_id] += 1

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
            increase_unique_id_occurence(unique_id)
            gene_ids = [str(val).lower() for col in gene_id_cols 
                        if pd.notna(val := row[col]) and val]
            for gene_id in gene_ids:  # Only add to dict if there are valid gene IDs
                pathways_dict[gene_id].append(f"{unique_id}_{pathway_occurences[unique_id]}")
    return dict(pathways_dict)

class TrieNode:
    def __init__(self):
        self.children = {}
        self.key = None  # Store the full key when a word ends here

class Trie:
    def __init__(self):
        self.root = TrieNode()

    def insert(self, key):
        """Insert a key into the trie."""
        node = self.root
        for char in key:
            if char not in node.children:
                node.children[char] = TrieNode()
            node = node.children[char]
        node.key = key  # Mark the end of the word and store the full key

    def search_substring(self, input_string):
        """
        Search if any substring of the input string matches a key in the trie.
        The substring must contain at least one English letter to be considered valid.
        
        Args:
            input_string: String to search within
        
        Returns:
            The matching key if found, None otherwise
        """
        def contains_letter(substr):
            """Check if substring contains at least one English letter."""
            return any(char in string.ascii_letters for char in substr)
        
        for start in range(len(input_string)):
            node = self.root
            current_substr = ""
            
            for char in input_string[start:]:
                current_substr += char
                if char in node.children:
                    node = node.children[char]
                    if node.key and contains_letter(current_substr):
                        return node.key
                else:
                    break
                    
        return None

# Process each file in a directory using multithreading
def process_file(file_path, pathway_dict, output_file, file_lock, trie):
    
    # Function to find the matching key from pathway_dict
    def find_key_in_id(id_value):
        id_lower = id_value.lower()
        match = trie.search_substring(id_lower)
        return match
    
    def get_pathway_from_dict(x, pathway_dict):
        key = find_key_in_id(x)
        if key in pathway_dict:
            return pathway_dict[key]
        return None
        
    print(f"Processing file: {file_path}")
    df = pd.read_csv(file_path)
    # num_genes = len(df)
    # Compile a regular expression pattern for all keys in pathway_dict
    # Apply the filtering and mapping
    new_df = df
    # Add source_file column
    new_df['source_file'] = file_path
    # Add pathway and index columns
    new_df['pathway']=new_df['id'].apply(lambda x: get_pathway_from_dict(x, pathway_dict))
    new_df['metabolic_gene'] = new_df['id'].apply(lambda x: find_key_in_id(x))

    new_df.to_csv(output_file, mode='w', header=True, index=False)
    print(f"Completed file: {file_path}")

# Process a whole genome directory
def process_genome_dir(genome_dir, pathway_dict, max_workers, output_dir, trie):
    file_paths = [os.path.join(genome_dir, f) for f in os.listdir(genome_dir) 
                  if f.endswith('.csv') and os.path.isfile(os.path.join(genome_dir, f))]
    
    output_paths = {f"{file_path}": os.path.join(output_dir, os.path.basename(file_path)) for file_path in file_paths}
    file_lock = Lock()
    with tqdm(total=len(file_paths), desc=f"Directory: {os.path.basename(genome_dir)}", unit="file") as pbar:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(
                    process_file, 
                    file_path, 
                    pathway_dict, 
                    output_paths[file_path],  
                    file_lock,
                    trie
                )
                for file_path in file_paths
            ]
            for future in as_completed(futures):
                pbar.update(1)
    print(f"Directory {genome_dir} completed")

def create_output_subdir(output_dir):
    """Create a subdirectory for each min_genes value."""
    subdir = os.path.join(output_dir, "metabolic_genes_with_letters")
    if not os.path.exists(subdir):
        os.makedirs(subdir)
    return subdir

def main():
    genome_dirs = [
        "/groups/itay_mayrose/alongonda/datasets/full_genomes/ensembl/processed_annotations_test_no_chloroplast_with_sequences",
        "/groups/itay_mayrose/alongonda/datasets/full_genomes/plaza/processed_annotations_with_chromosomes_no_chloroplast_with_sequences",
        "/groups/itay_mayrose/alongonda/datasets/full_genomes/phytozome/processed_annotations_with_chromosomes_no_chloroplast_with_sequences"
    ]
    pathways_file = "/groups/itay_mayrose/alongonda/datasets/plantcyc/all_organisms/merged_pathways.csv"
        
    max_workers = min(32, os.cpu_count())
    print(f"Using {max_workers} workers")
    print("Loading pathways...")

    # Load pathway dictionary and build Aho-Corasick automaton
    pathway_dict = create_pathways_dict(pathways_file)
    
    # Build the trie with pathway_dict keys
    trie = Trie()
    for key in pathway_dict.keys():
        trie.insert(key.lower())

    for genome_dir in genome_dirs:
        print(f"Processing genome directory: {genome_dir}")
        output_dir = create_output_subdir(genome_dir)
        process_genome_dir(genome_dir, pathway_dict, max_workers, output_dir, trie)
            
if __name__ == "__main__":
    main()