#!/usr/bin/env python3
"""
Random MGC Sampler
Creates 10,000 random gene groups from annotated genome files
Each group contains 3-5 randomly selected genes from a single organism
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path
import random
from tqdm import tqdm
import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from multiprocessing import cpu_count
import time

def get_annotated_files(annotated_dir):
    """Get all annotated CSV files"""
    annotated_path = Path(annotated_dir)
    csv_files = list(annotated_path.glob("*.csv"))
    print(f"Found {len(csv_files)} annotated files")
    return csv_files

def sample_genes_from_file(csv_file, group_size, sample_index):
    """Sample genes from a single annotated file"""
    try:
        df = pd.read_csv(csv_file)
        
        # Check required columns
        required_cols = ['id', 'kegg_ids', 'sequence']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"Warning: Missing columns {missing_cols} in {csv_file.name}")
            return None
        
        # Filter rows with valid data
        valid_df = df.dropna(subset=['id', 'sequence'])
        valid_df = valid_df[valid_df['sequence'].str.len() > 20]  # Filter short sequences
        
        if len(valid_df) < group_size:
            print(f"Warning: {csv_file.name} has only {len(valid_df)} valid genes, need {group_size}")
            return None
        
        # Sample genes
        sampled = valid_df.sample(n=group_size, random_state=None)
        
        # Create result format with RANDOM_MGC_X basename
        result = []
        basename = f"RANDOM_MGC_{sample_index}"
        
        for _, row in sampled.iterrows():
            result.append({
                'mgc_candidate': basename,
                'gene_id': row['id'],
                'kegg_id': row.get('kegg_ids', ''),
                'source_file': str(csv_file.absolute()),
                'sequence': row['sequence']
            })
        
        return result
        
    except Exception as e:
        print(f"Error processing {csv_file.name}: {e}")
        return None

def process_file_batch(file_batch_info):
    """Process a batch of samples for concurrent execution"""
    csv_file, num_samples_for_file, start_index = file_batch_info
    results = []
    
    for i in range(num_samples_for_file):
        sample_index = start_index + i
        group_size = random.choice([3, 4, 5])
        
        sampled_genes = sample_genes_from_file(csv_file, group_size, sample_index)
        
        if sampled_genes:
            results.append((sample_index, sampled_genes))
    
    return results

def create_random_mgc_samples(annotated_dir, output_dir, target_samples=10000, max_workers=None):
    """Create 10,000 random MGC samples"""
    
    print(f"Creating {target_samples} random MGC samples...")
    print(f"Input directory: {annotated_dir}")
    print(f"Output directory: {output_dir}")
    
    # Setup
    os.makedirs(output_dir, exist_ok=True)
    random.seed(42)  # For reproducibility
    np.random.seed(42)
    
    # Get all annotated files
    csv_files = get_annotated_files(annotated_dir)
    if not csv_files:
        print("No annotated CSV files found!")
        return
    
    # Calculate how many samples to generate from each file
    samples_per_file = target_samples // len(csv_files)
    extra_samples = target_samples % len(csv_files)
    
    print(f"Generating ~{samples_per_file} samples per file")
    if extra_samples > 0:
        print(f"Extra {extra_samples} samples from random files")
    
    # Determine number of workers
    if max_workers is None:
        max_workers = min(cpu_count(), len(csv_files), 8)  # Cap at 8 to avoid overwhelming I/O
    
    print(f"Using {max_workers} worker processes for parallel processing")
    
    # Prepare batch information for concurrent processing
    file_batches = []
    current_index = 1
    
    for i, csv_file in enumerate(csv_files):
        # Determine number of samples for this file
        num_samples = samples_per_file
        if i < extra_samples:  # Distribute extra samples
            num_samples += 1
        
        if num_samples > 0:
            file_batches.append((csv_file, num_samples, current_index))
            current_index += num_samples
    
    # Process files concurrently
    start_time = time.time()
    successful_samples = 0
    total_samples_attempted = 0
    
    print("Starting concurrent processing...")
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all batch jobs
        future_to_file = {executor.submit(process_file_batch, batch): batch[0] 
                         for batch in file_batches}
        
        # Process results as they complete
        for future in tqdm(as_completed(future_to_file), total=len(file_batches), desc="Processing files"):
            csv_file = future_to_file[future]
            try:
                batch_results = future.result()
                total_samples_attempted += len([r for r in batch_results if r])
                
                # Save each successful sample
                for sample_index, sampled_genes in batch_results:
                    if sampled_genes:
                        output_file = os.path.join(output_dir, f"RANDOM_MGC_{sample_index}.csv")
                        sample_df = pd.DataFrame(sampled_genes)
                        sample_df.to_csv(output_file, index=False)
                        successful_samples += 1
                
                # Progress update
                if successful_samples % 1000 == 0 and successful_samples > 0:
                    elapsed = time.time() - start_time
                    rate = successful_samples / elapsed if elapsed > 0 else 0
                    print(f"Created {successful_samples} samples ({rate:.1f} samples/sec)")
                    
            except Exception as exc:
                print(f'File {csv_file.name} generated an exception: {exc}')
    
    # Final timing and statistics
    total_time = time.time() - start_time
    
    print(f"\n✅ Random MGC sampling complete!")
    print(f"Total samples attempted: {total_samples_attempted}")
    print(f"Successful samples: {successful_samples}")
    print(f"Success rate: {successful_samples/total_samples_attempted*100:.1f}%" if total_samples_attempted > 0 else "Success rate: 0%")
    print(f"Processing time: {total_time:.1f} seconds")
    print(f"Average rate: {successful_samples/total_time:.1f} samples/second" if total_time > 0 else "")
    print(f"Output directory: {output_dir}")

def main():
    parser = argparse.ArgumentParser(description='Create random MGC samples from annotated genomes')
    parser.add_argument('--annotated_dir', 
                       default='/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/annotated_genomes_metabolic',
                       help='Directory containing annotated CSV files')
    parser.add_argument('--output_dir',
                       default='/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/random_mgc_candidates_csv_files',
                       help='Output directory for random samples')
    parser.add_argument('--num_samples', type=int, default=10000,
                       help='Number of random samples to create')
    parser.add_argument('--max_workers', type=int, default=None,
                       help='Maximum number of worker processes (default: auto-detect)')
    
    args = parser.parse_args()
    
    # Validate input directory
    if not os.path.exists(args.annotated_dir):
        print(f"Error: Annotated directory not found: {args.annotated_dir}")
        return
    
    # Create random samples
    create_random_mgc_samples(args.annotated_dir, args.output_dir, args.num_samples, args.max_workers)

if __name__ == "__main__":
    main()