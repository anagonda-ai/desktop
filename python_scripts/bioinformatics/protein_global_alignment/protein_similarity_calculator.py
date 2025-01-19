from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import numpy as np
import argparse
import os
import sys
from itertools import combinations
import pandas as pd
import datetime
import concurrent.futures
from tqdm import tqdm
import multiprocessing
import gc
from queue import Queue
from threading import Thread
import math

class SequenceCache:
    """
    Cache for protein sequences to avoid repeated file reads
    """
    def __init__(self, max_size=10):
        self.cache = {}
        self.max_size = max_size
        
    def get_sequences(self, filename):
        if filename not in self.cache:
            if len(self.cache) >= self.max_size:
                # Remove least recently used item
                self.cache.pop(next(iter(self.cache)))
            self.cache[filename] = read_protein_sequences(filename)
        return self.cache[filename]

def read_protein_sequences(filename):
    """
    Read protein sequences from a FASTA file
    """
    if not os.path.exists(filename):
        sys.exit(f"Error: File {filename} does not exist")
        
    sequences = {}
    for record in SeqIO.parse(filename, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def calculate_global_alignment_score(seq1, seq2):
    """
    Calculate global alignment score between two protein sequences
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    
    actual_score = aligner.score(seq1, seq2)
    min_blosum_score = -4
    worst_score = min_blosum_score * max(len(seq1), len(seq2))
    longer_seq = seq1 if len(seq1) >= len(seq2) else seq2
    best_score = aligner.score(longer_seq, longer_seq)
    
    normalized_score = (actual_score - worst_score) / (best_score - worst_score)
    return max(0.0, min(1.0, normalized_score))

def process_pair(args):
    """
    Process a single sequence pair
    """
    seq1_id, seq1, seq2_id, seq2, file1, file2 = args
    try:
        score = calculate_global_alignment_score(seq1, seq2)
        return {
            'File1': file1,
            'File2': file2,
            'Sequence1': seq1_id,
            'Sequence2': seq2_id,
            'Similarity_Score': score
        }
    except Exception as e:
        return None

def process_mega_chunk(chunk_data, sequence_cache):
    """
    Process a large chunk of sequence pairs with internal parallelization
    """
    file1_path, file2_path, start_idx, end_idx, inner_chunk_size = chunk_data
    
    try:
        # Get sequences from cache
        group1_sequences = sequence_cache.get_sequences(file1_path)
        group2_sequences = sequence_cache.get_sequences(file2_path)
        
        file1 = os.path.basename(file1_path)
        file2 = os.path.basename(file2_path)
        
        # Prepare sequence pairs
        pairs = []
        pairs_done = 0
        
        for seq1_id, seq1 in group1_sequences.items():
            for seq2_id, seq2 in group2_sequences.items():
                if start_idx <= pairs_done < end_idx:
                    pairs.append((seq1_id, seq1, seq2_id, seq2, file1, file2))
                pairs_done += 1
                if pairs_done >= end_idx:
                    break
            if pairs_done >= end_idx:
                break
        
        # Process pairs in parallel using ThreadPoolExecutor
        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
            futures = [executor.submit(process_pair, pair) for pair in pairs]
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                if result:
                    results.append(result)
        
        return results
    
    except Exception as e:
        print(f"Error processing mega chunk: {str(e)}")
        return []

def parallel_write_results(queue, output_file):
    """
    Continuously write results from queue to file
    """
    while True:
        results = queue.get()
        if results is None:  # Poison pill
            break
        
        df = pd.DataFrame(results)
        df.to_csv(output_file, mode='a', header=False, index=False)
        del df
        gc.collect()

def perform_all_vs_all_alignment(input_dir, output_dir, max_workers=None, mega_chunk_size=10000):
    """
    Perform all-vs-all alignments with enhanced parallelization
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if max_workers is None:
        max_workers = max(1, multiprocessing.cpu_count() - 2)
    
    # Initialize sequence cache
    sequence_cache = SequenceCache(max_size=20)
    
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    detailed_output = os.path.join(output_dir, f'detailed_alignments_{timestamp}.csv')
    summary_output = os.path.join(output_dir, f'similarity_matrix_{timestamp}.csv')
    
    # Initialize output file with headers
    pd.DataFrame(columns=['File1', 'File2', 'Sequence1', 'Sequence2', 'Similarity_Score'])\
        .to_csv(detailed_output, index=False)
    
    # Set up result writing queue and thread
    result_queue = Queue(maxsize=100)
    writer_thread = Thread(target=parallel_write_results, args=(result_queue, detailed_output))
    writer_thread.start()
    
    print(f"Processing files using {max_workers} workers...")
    
    # Get file pairs
    fasta_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.fasta')]
    file_pairs = list(combinations(fasta_files, 2))
    
    # Process file pairs with parallel mega chunks
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        for pair_idx, (file1_path, file2_path) in enumerate(file_pairs, 1):
            try:
                # Calculate total pairs
                seq_count1 = len(list(SeqIO.parse(file1_path, "fasta")))
                seq_count2 = len(list(SeqIO.parse(file2_path, "fasta")))
                total_pairs = seq_count1 * seq_count2
                
                print(f"\nProcessing file pair {pair_idx}/{len(file_pairs)}: "
                      f"{os.path.basename(file1_path)} vs {os.path.basename(file2_path)}")
                
                # Create mega chunks
                mega_chunks = []
                for start in range(0, total_pairs, mega_chunk_size):
                    end = min(start + mega_chunk_size, total_pairs)
                    mega_chunks.append((file1_path, file2_path, start, end, 1000))
                
                # Process mega chunks in parallel
                futures = [executor.submit(process_mega_chunk, chunk, sequence_cache) 
                          for chunk in mega_chunks]
                
                for future in tqdm(concurrent.futures.as_completed(futures), 
                                 total=len(futures),
                                 desc="Processing mega chunks"):
                    results = future.result()
                    if results:
                        result_queue.put(results)
                
            except Exception as e:
                print(f"Error processing file pair: {str(e)}")
                continue
    
    # Signal writer thread to finish
    result_queue.put(None)
    writer_thread.join()
    
    # Create summary matrix
    print("\nGenerating summary matrix...")
    try:
        df = pd.read_csv(detailed_output)
        summary_df = df.groupby(['File1', 'File2'])['Similarity_Score'].mean().unstack()
        summary_df.to_csv(summary_output)
    except Exception as e:
        print(f"Error generating summary matrix: {str(e)}")
    
    return detailed_output, summary_output

def main():
    parser = argparse.ArgumentParser(description='Perform all-vs-all protein similarity analysis')
    parser.add_argument('-i', '--input_dir', required=True, help='Directory containing FASTA files')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory for output files')
    parser.add_argument('-w', '--workers', type=int, default=None, 
                        help='Number of worker processes (default: CPU cores - 2)')
    parser.add_argument('-c', '--chunk_size', type=int, default=10000,
                        help='Size of mega chunks (default: 10000)')
    
    args = parser.parse_args()
    
    try:
        detailed_file, summary_file = perform_all_vs_all_alignment(
            args.input_dir,
            args.output_dir,
            args.workers,
            args.chunk_size
        )
        print(f"\nAnalysis complete!")
        print(f"Detailed results saved to: {detailed_file}")
        print(f"Summary matrix saved to: {summary_file}")
                
    except Exception as e:
        sys.exit(f"Error: {str(e)}")

if __name__ == "__main__":
    main()