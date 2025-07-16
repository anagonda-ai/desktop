#!/bin/bash
#SBATCH --job-name=kegg_annot
#SBATCH --output=/groups/itay_mayrose/alongonda/desktop/slurm_logs/kegg_%j.out
#SBATCH --error=/groups/itay_mayrose/alongonda/desktop/slurm_logs/kegg_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --partition=itaym
#SBATCH --time=3-24:00:00

annotated_file=$1
output_file=$2
window_size=$3
min_genes=$4

python3 - <<EOF
import pandas as pd
import numpy as np
import os
from itertools import combinations
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
from collections import defaultdict, Counter
import gc

annotated_file = "$annotated_file"
output_file = "$output_file"
window_size = int("$window_size")
min_genes = int("$min_genes")
file_lock = threading.Lock()

class PathwayMatcher:
    def __init__(self, min_genes):
        self.min_genes = min_genes
        self.pathway_cache = {}
        
    def find_common_pathway_vectorized(self, gene_ids, pathway_arrays):
        """Ultra-fast pathway finding using numpy operations and caching"""
        # Create cache key
        cache_key = tuple(sorted(gene_ids))
        if cache_key in self.pathway_cache:
            return self.pathway_cache[cache_key]
        
        if len(gene_ids) < self.min_genes:
            result = (None, [])
            self.pathway_cache[cache_key] = result
            return result
        
        # Get pathway sets for all genes
        valid_genes = []
        pathway_sets = []
        
        for gene_id in gene_ids:
            if gene_id in pathway_arrays:
                pathways = pathway_arrays[gene_id]
                if pathways:
                    valid_genes.append(gene_id)
                    pathway_sets.append(pathways)
        
        if len(valid_genes) < self.min_genes:
            result = (None, [])
            self.pathway_cache[cache_key] = result
            return result
        
        # Find intersection efficiently using Counter
        pathway_counts = Counter()
        for pathway_set in pathway_sets:
            for pathway in pathway_set:
                pathway_counts[pathway] += 1
        
        # Find pathways present in at least min_genes
        common_pathways = [pathway for pathway, count in pathway_counts.items() 
                          if count >= self.min_genes]
        
        if not common_pathways:
            result = (None, [])
            self.pathway_cache[cache_key] = result
            return result
        
        # Return first common pathway and genes that have it
        target_pathway = common_pathways[0]
        genes_with_pathway = [gene for gene in valid_genes 
                             if target_pathway in pathway_arrays[gene]][:self.min_genes]
        
        result = (target_pathway, genes_with_pathway)
        self.pathway_cache[cache_key] = result
        return result

def preprocess_annotations(df):
    """Vectorized preprocessing of pathway and annotation data"""
    # Convert pathways to sets in vectorized manner
    pathway_dict = {}
    annotation_dict = {}
    
    # Process in chunks to avoid memory issues
    chunk_size = 10000
    for i in range(0, len(df), chunk_size):
        chunk = df.iloc[i:i+chunk_size]
        
        for _, row in chunk.iterrows():
            gene_id = row['id']
            
            # Process pathways
            if pd.notna(row['pathway']):
                pathways = frozenset(p.strip() for p in row['pathway'].split(",") if p.strip())
                pathway_dict[gene_id] = pathways
            else:
                pathway_dict[gene_id] = frozenset()
            
            # Process annotations
            if pd.notna(row['annotation']):
                annotations = [a.strip() for a in row['annotation'].split(",") if a.strip()]
                annotation_dict[gene_id] = annotations[0] if annotations else ''
            else:
                annotation_dict[gene_id] = ''
    
    return pathway_dict, annotation_dict

def process_chromosome_vectorized(chr_data, window_size, min_genes, annotated_file, pathway_matcher):
    """Highly optimized chromosome processing with numpy vectorization"""
    if len(chr_data) < min_genes:
        return []
    
    # Convert to numpy arrays for maximum speed
    positions = chr_data['index'].values
    gene_ids = chr_data['id'].values
    starts = chr_data['start'].values
    ends = chr_data['end'].values
    
    # Preprocess annotations once
    pathway_dict, annotation_dict = preprocess_annotations(chr_data)
    
    matches = []
    processed_combinations = set()
    last_match = None
    
    # Use sliding window with numpy operations
    n_genes = len(positions)
    
    # Pre-compute all possible windows using numpy
    for i in range(n_genes):
        start_pos = positions[i]
        
        # Find window boundaries using vectorized operations
        within_window = positions >= start_pos
        within_window &= positions <= start_pos + window_size
        window_indices = np.where(within_window)[0]
        
        # Skip if window too small
        if len(window_indices) < min_genes:
            continue
        
        # Get genes in current window
        window_genes = gene_ids[window_indices]
        
        # Find common pathway
        pathway, metabolic_genes = pathway_matcher.find_common_pathway_vectorized(
            window_genes, pathway_dict
        )
        
        if not pathway or not metabolic_genes:
            continue
        
        # Check if this combination was already processed
        combo_key = tuple(sorted(metabolic_genes))
        if combo_key in processed_combinations:
            continue
        processed_combinations.add(combo_key)
        
        # Get annotations for metabolic genes
        metabolic_annotations = [annotation_dict.get(gene, '') for gene in metabolic_genes]
        
        # Calculate window boundaries
        window_start = starts[window_indices].min()
        window_end = ends[window_indices].max()
        
        current_match = {
            'pathway': pathway,
            'genes': window_genes.tolist(),
            'metabolic_genes': metabolic_genes,
            'metabolic_genes_annotations': metabolic_annotations,
            'start': window_start,
            'end': window_end
        }
        
        # Optimize merging logic
        should_merge = (
            last_match and 
            pathway == last_match['pathway'] and 
            len(metabolic_genes) >= 2 and 
            len(last_match['metabolic_genes']) >= 2 and
            metabolic_genes[:2] == last_match['metabolic_genes'][-2:]
        )
        
        if should_merge:
            # Use set operations for faster merging
            last_genes_set = set(last_match['genes'])
            last_metabolic_set = set(last_match['metabolic_genes'])
            
            new_genes = [g for g in current_match['genes'] if g not in last_genes_set]
            new_metabolic_genes = [g for g in metabolic_genes if g not in last_metabolic_set]
            new_annotations = [ann for gene, ann in zip(metabolic_genes, metabolic_annotations)
                             if gene in new_metabolic_genes]
            
            last_match['genes'].extend(new_genes)
            last_match['metabolic_genes'].extend(new_metabolic_genes)
            last_match['metabolic_genes_annotations'].extend(new_annotations)
            last_match['end'] = max(last_match['end'], current_match['end'])
        else:
            if last_match:
                matches.append(format_match(last_match, annotated_file))
            last_match = current_match
    
    # Add final match
    if last_match:
        matches.append(format_match(last_match, annotated_file))
    
    return matches

def format_match(match, annotated_file):
    """Format match for output"""
    return {
        'pathway': match['pathway'],
        'genes': ','.join(match['genes']),
        'metabolic_genes': ','.join(match['metabolic_genes']),
        'metabolic_genes_annotations': ','.join(match['metabolic_genes_annotations']),
        'start': match['start'],
        'end': match['end'],
        'source_file': annotated_file
    }

def process_annotated_file_optimized(annotated_file, output_file, window_size, min_genes):
    """Main processing with maximum optimization"""
    print(f"Processing: {annotated_file}")
    
    try:
        # Read with optimized settings
        df = pd.read_csv(
            annotated_file, 
            dtype={'chromosome': 'category'},
            engine='c',  # Use C parser for speed
            memory_map=True  # Use memory mapping for large files
        )
    except Exception as e:
        print(f"Error reading {annotated_file}: {e}")
        return
    
    # Add index and filter in one operation
    df.reset_index(drop=True, inplace=True)
    df["index"] = df.index
    
    # Filter more efficiently
    valid_mask = df['pathway'].notna() & (df['pathway'].str.len() > 0)
    filtered_df = df[valid_mask].copy()
    
    if filtered_df.empty:
        print(f"⚠️ {annotated_file} - No valid pathways found")
        return
    
    # Create pathway matcher instance
    pathway_matcher = PathwayMatcher(min_genes)
    
    # Group chromosomes and optimize memory usage
    chromosome_groups = list(filtered_df.groupby('chromosome', observed=True))
    
    # Process with optimal thread count
    max_workers = min(16, len(chromosome_groups), os.cpu_count())
    all_matches = []
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit tasks with progress tracking
        future_to_chr = {}
        for chr_name, chr_data in chromosome_groups:
            # Sort by position for optimal sliding window
            chr_data_sorted = chr_data.sort_values('index').reset_index(drop=True)
            future = executor.submit(
                process_chromosome_vectorized, 
                chr_data_sorted, window_size, min_genes, annotated_file, pathway_matcher
            )
            future_to_chr[future] = chr_name
        
        # Process results with progress indication
        completed = 0
        total = len(future_to_chr)
        
        for future in as_completed(future_to_chr):
            chr_name = future_to_chr[future]
            completed += 1
            
            try:
                matches = future.result()
                all_matches.extend(matches)
                if matches:
                    print(f"  Chr {chr_name}: {len(matches)} matches ({completed}/{total})")
            except Exception as e:
                print(f"  Error processing chromosome {chr_name}: {e}")
    
    # Efficient output writing
    if all_matches:
        write_results_optimized(all_matches, output_file)
        print(f"✔️ {annotated_file} - Total matches: {len(all_matches)}")
    else:
        print(f"⚠️ {annotated_file} - No matches found")
    
    # Force garbage collection
    del df, filtered_df, all_matches
    gc.collect()

def write_results_optimized(matches, output_file):
    """Optimized result writing with batching"""
    with file_lock:
        file_exists = os.path.exists(output_file)
        
        # Write in chunks to avoid memory issues with large result sets
        chunk_size = 10000
        for i in range(0, len(matches), chunk_size):
            chunk = matches[i:i+chunk_size]
            chunk_df = pd.DataFrame(chunk)
            
            # Write header only for first chunk of first file
            write_header = not file_exists and i == 0
            chunk_df.to_csv(
                output_file, 
                mode='a' if file_exists or i > 0 else 'w',
                header=write_header,
                index=False,
                chunksize=1000  # Write in smaller chunks
            )
            file_exists = True  # After first write, file exists

# Execute main function
process_annotated_file_optimized(annotated_file, output_file, window_size, min_genes)
EOF