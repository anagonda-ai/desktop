import os
import csv
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from threading import Lock, Event
import hashlib

def create_pathways_dict(pathways_file):
    """Read pathway file using chunks to avoid memory issues"""
    pathways_dict = defaultdict(set)
    
    chunk_size = 1000
    
    for chunk in pd.read_csv(pathways_file, chunksize=chunk_size, low_memory=False):
        gene_id_cols = [col for col in chunk.columns if col.startswith('GENE-ID.')]
        required_cols = ['UNIQUE-ID'] + gene_id_cols
        
        chunk = chunk[required_cols]
        
        for _, row in chunk.iterrows():
            unique_id = row['UNIQUE-ID']
            gene_ids = {str(val).lower() for col in gene_id_cols 
                        if pd.notna(val := row[col]) and val}
            if gene_ids:
                pathways_dict[unique_id].update(gene_ids)
    
    return dict(pathways_dict)

def read_csv_files(input_csv_dir):
    return [os.path.join(input_csv_dir, f) for f in os.listdir(input_csv_dir)
            if f.endswith('.csv') and os.path.isfile(os.path.join(input_csv_dir, f))]

def generate_group_hash(group):
    """
    Generate a unique hash for a group to identify duplicates
    Sorts genes to ensure consistent hashing
    """
    # Sort genes to make hash consistent regardless of order
    sorted_genes = sorted(group['cluster_genes'].split(','))
    
    # Create a hash string that captures key identifying information
    hash_string = f"{group['pathway']}|{'|'.join(sorted_genes)}"
    
    # Use MD5 to create a fixed-length unique identifier
    return hashlib.md5(hash_string.encode()).hexdigest()

class UniqueMatchTracker:
    def __init__(self):
        self._seen_matches = set()
        self._lock = Lock()
    
    def is_unique_match(self, pathway_cluster):
        """
        Check if the match is unique and add it to seen matches if so
        
        Args:
        - group: Dictionary representing a pathway and a cluster of genes
        
        Returns:
        - Boolean indicating if this is a new, unique match
        """
        group_hash = generate_group_hash(pathway_cluster)
        
        with self._lock:
            if group_hash in self._seen_matches:
                return False
            self._seen_matches.add(group_hash)
            return True

def process_window(window_data, pathway_dict, output_file, file_lock, 
                   file_path, unique_tracker, min_genes=3):
    """
    Process a window of genes and write unique matches
    """
    gene_ids = window_data['id'].tolist()
    matches_found = 0
    
    for pathway, pathway_genes in pathway_dict.items():
        # Find matches between window genes and pathway genes
        matches = [g for g in pathway_genes if any(g.lower() in gene_id.lower() for gene_id in gene_ids)]
        
        # Check if we have enough matches
        if matches and len(matches) >= min_genes:
            # Prepare group information
            group = {
                'pathway': pathway,
                'genes': ','.join(gene_ids),
                'cluster_genes': ','.join(sorted(matches)),
                'start': window_data['start'].min(),
                'end': window_data['end'].max(),
                'source_file': os.path.basename(file_path)
            }
            
            pathway_cluster = {
                'pathway': pathway,
                'cluster_genes': ','.join(sorted(matches))}
            
            # Check if this is a unique match
            if unique_tracker.is_unique_match(pathway_cluster):
                matches_found += 1
                
                # Thread-safe file writing
                with file_lock:
                    mode = 'a' if os.path.exists(output_file) else 'w'
                    header = not os.path.exists(output_file)
                    
                    df_group = pd.DataFrame([group])
                    df_group.to_csv(output_file, mode=mode, header=header, index=False)
                
                # Print match details
                print(f"Unique Match Found: Pathway={pathway}, " + 
                      f"Genes={','.join(gene_ids)}, " + 
                      f"Matching Genes={','.join(matches)}, " + 
                      f"File={os.path.basename(file_path)}")
    
    return matches_found

def process_file(file_path, pathway_dict, output_file, file_lock, 
                 unique_tracker, window_size=10, min_genes=3, chunk_size=1000):
    print(f"Processing file: {file_path}")
    total_matches = 0
    
    try:
        # Use chunksize to read file incrementally
        csv_reader = pd.read_csv(file_path, chunksize=chunk_size, usecols=['id', 'start', 'end'])
        
        for chunk_idx, chunk in enumerate(csv_reader):
            # Sort chunk by start and end
            chunk = chunk.sort_values(by=['start', 'end'])
            
            # Process windows within the chunk
            for j in range(len(chunk) - window_size + 1):
                window = chunk.iloc[j:j + window_size]
                chunk_matches = process_window(
                    window_data=window, 
                    pathway_dict=pathway_dict, 
                    output_file=output_file, 
                    file_lock=file_lock, 
                    file_path=file_path,
                    unique_tracker=unique_tracker,
                    min_genes=min_genes
                )
                total_matches += chunk_matches
            
            # Print progress for every 10 chunks
            if chunk_idx % 10 == 0:
                print(f"Processed chunk {chunk_idx} in {os.path.basename(file_path)}")
        
        print(f"File {os.path.basename(file_path)} - Total Unique Matches: {total_matches}")
                
    except Exception as e:
        print(f"Error processing {file_path}: {str(e)}")
        return 0
    
    return total_matches

def process_genome_dir(genome_dir, pathway_dict, output_file, max_workers):
    file_paths = read_csv_files(genome_dir)
    
    # Use locks to prevent race conditions
    file_lock = Lock()
    unique_tracker = UniqueMatchTracker()
    
    total_dir_matches = 0
    
    # Process files in smaller batches to manage memory
    batch_size = 4
    for i in range(0, len(file_paths), batch_size):
        batch_paths = file_paths[i:i + batch_size]
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Create futures for batch processing
            futures = [
                executor.submit(
                    process_file, 
                    fp, 
                    pathway_dict, 
                    output_file, 
                    file_lock,
                    unique_tracker
                ) 
                for fp in batch_paths
            ]
            
            # Wait for all futures to complete and collect matches
            for future in as_completed(futures):
                try:
                    batch_matches = future.result()
                    total_dir_matches += batch_matches
                except Exception as exc:
                    print(f"Unexpected error in processing: {exc}")
        
        print(f"Completed batch {i//batch_size + 1}/{(len(file_paths) + batch_size - 1)//batch_size}")
    
    print(f"Directory {os.path.basename(genome_dir)} - Total Unique Matches: {total_dir_matches}")
    return total_dir_matches

def main():
    genome_dirs = [
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/ensembl/processed_annotations_sorted",
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/plaza/processed_annotations_sorted",
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/phytozome/processed_annotations"
    ]
    pathways_file = "/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms/merged_pathways.csv"
    output_file = "/groups/itay_mayrose/alongonda/desktop/potential_groups_immediate.csv"

    # Calculate optimal number of workers
    max_workers = min(32, os.cpu_count())
    print(f"Using {max_workers} workers")
    
    print("Loading pathway dictionary...")
    pathway_dict = create_pathways_dict(pathways_file)
    print(f"Loaded {len(pathway_dict)} pathways")
    
    # Remove output file if it exists to start fresh
    if os.path.exists(output_file):
        os.remove(output_file)
    
    total_global_matches = 0
    for genome_dir in genome_dirs:
        print(f"Processing directory: {genome_dir}")
        dir_matches = process_genome_dir(genome_dir, pathway_dict, output_file, max_workers)
        total_global_matches += dir_matches
    
    print(f"TOTAL UNIQUE MATCHES FOUND: {total_global_matches}")
    print(f"Results saved to: {output_file}")

if __name__ == "__main__":
    main()