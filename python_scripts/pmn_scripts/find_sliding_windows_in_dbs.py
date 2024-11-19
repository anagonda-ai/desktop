import os
import csv
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from itertools import islice

def create_pathways_dict(pathways_file):
    """Read pathway file using chunks to avoid memory issues"""
    pathways_dict = defaultdict(set)
    
    # Read CSV in chunks with low_memory=False
    chunk_size = 1000
    
    for chunk in pd.read_csv(pathways_file, 
                             chunksize=chunk_size, 
                             low_memory=False):
        # Dynamically select columns with prefix "GENE-ID."
        gene_id_cols = [col for col in chunk.columns if col.startswith('GENE-ID.')]
        required_cols = ['UNIQUE-ID'] + gene_id_cols
        
        # Filter the chunk to include only the required columns
        chunk = chunk[required_cols]
        
        for _, row in chunk.iterrows():
            unique_id = row['UNIQUE-ID']
            # Get all gene IDs, filtering out NaN and empty values
            gene_ids = {str(val).lower() for col in gene_id_cols 
                        if pd.notna(val := row[col]) and val}
            if gene_ids:
                pathways_dict[unique_id].update(gene_ids)
    
    return dict(pathways_dict)

def read_csv_files(input_csv_dir):
    return [os.path.join(input_csv_dir, f) for f in os.listdir(input_csv_dir)
            if f.endswith('.csv') and os.path.isfile(os.path.join(input_csv_dir, f))]

def process_window(window_data, pathway_dict, min_genes=3):
    print(f"Processing window: {window_data['start'].min()} - {window_data['end'].max()}")
    potential_groups = set()
    gene_ids = window_data['id'].tolist()
    
    # Pre-compute gene matches for all pathways
    gene_matches = defaultdict(set)
    for gene_id in gene_ids:
        gene_id_lower = gene_id.lower()
        for pathway, pathway_genes in pathway_dict.items():
            matches = {g for g in pathway_genes if g.lower() in gene_id_lower}
            if matches:
                gene_matches[pathway].update(matches)
    
    # Only process pathways with sufficient matches
    for pathway, matches in gene_matches.items():
        if len(matches) >= min_genes:
            group = {
                'pathway': pathway,
                'genes': tuple(sorted(gene_ids)),
                'cluster_genes': tuple(sorted(matches)),
                'start': window_data['start'].min(),
                'end': window_data['end'].max()
            }
            potential_groups.add(tuple(sorted(group.items())))
    
    return potential_groups

def process_file(file_path, pathway_dict, window_size=10, min_genes=3):
    print(f"Processing file: {file_path}")
    potential_groups = set()
    
    try:
        # Read only necessary columns
        df = pd.read_csv(file_path, usecols=['id', 'start', 'end'])
        
        # Process in smaller chunks to manage memory
        chunk_size = window_size * 100
        for i in range(0, len(df), chunk_size):
            chunk = df.iloc[i:i + chunk_size].sort_values(by=['start', 'end'])
            
            # Process windows within the chunk
            for j in range(len(chunk) - window_size + 1):
                window = chunk.iloc[j:j + window_size]
                window_groups = process_window(window, pathway_dict, min_genes)
                potential_groups.update(window_groups)
                
            if i % 1000 == 0:
                print(f"Processed {i}/{len(df)} rows in {file_path}")
                
    except Exception as e:
        print(f"Error processing {file_path}: {str(e)}")
        return set()
    
    return potential_groups

def save_potential_groups_to_file(potential_groups, output_file):
    # Save in chunks to handle large datasets
    chunk_size = 1000
    rows = []
    
    for group in potential_groups:
        group_dict = dict(group)
        rows.append({
            'pathway': group_dict['pathway'],
            'genes': ','.join(group_dict['genes']),
            'cluster_genes': ','.join(group_dict['cluster_genes']),
            'start': group_dict['start'],
            'end': group_dict['end']
        })
        
        if len(rows) >= chunk_size:
            mode = 'w' if not os.path.exists(output_file) else 'a'
            header = not os.path.exists(output_file)
            pd.DataFrame(rows).to_csv(output_file, mode=mode, header=header, index=False)
            rows = []
    
    # Save any remaining rows
    if rows:
        mode = 'w' if not os.path.exists(output_file) else 'a'
        header = not os.path.exists(output_file)
        pd.DataFrame(rows).to_csv(output_file, mode=mode, header=header, index=False)

def process_genome_dir(genome_dir, pathway_dict, max_workers):
    file_paths = read_csv_files(genome_dir)
    all_potential_groups = set()
    
    # Process files in smaller batches to manage memory
    batch_size = 4
    for i in range(0, len(file_paths), batch_size):
        batch_paths = file_paths[i:i + batch_size]
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_file = {executor.submit(process_file, fp, pathway_dict): fp 
                            for fp in batch_paths}
            
            for future in as_completed(future_to_file):
                try:
                    potential_groups = future.result()
                    all_potential_groups.update(potential_groups)
                except Exception as exc:
                    file_path = future_to_file[future]
                    print(f"Error processing {file_path}: {exc}")
        
        print(f"Completed batch {i//batch_size + 1}/{(len(file_paths) + batch_size - 1)//batch_size}")
    
    return all_potential_groups

def main():
    genome_dirs = [
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/ensembl/processed_annotations_sorted",
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/plaza/processed_annotations_sorted",
        "/groups/itay_mayrose_nosnap/alongonda/full_genomes/phytozome/processed_annotations"
    ]
    pathways_file = "/groups/itay_mayrose/alongonda/desktop/plantcyc/all_organisms/merged_pathways.csv"
    output_file = "/groups/itay_mayrose/alongonda/desktop/potential_groups.csv"

    # Calculate optimal number of workers
    max_workers = min(32, os.cpu_count())  # Reduced from 32 to prevent memory issues
    print(f"Using {max_workers} workers")
    
    print("Loading pathway dictionary...")
    pathway_dict = create_pathways_dict(pathways_file)
    print(f"Loaded {len(pathway_dict)} pathways")
    
    all_potential_groups = set()
    for genome_dir in genome_dirs:
        print(f"Processing directory: {genome_dir}")
        dir_groups = process_genome_dir(genome_dir, pathway_dict, max_workers)
        all_potential_groups.update(dir_groups)
        print(f"Found {len(dir_groups)} groups in {genome_dir}")
    
    print(f"Saving {len(all_potential_groups)} potential groups...")
    save_potential_groups_to_file(all_potential_groups, output_file)
    print(f"Results saved to: {output_file}")

if __name__ == "__main__":
    main()