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
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import multiprocessing
from queue import Queue
from threading import Thread
from functools import lru_cache
import warnings
warnings.filterwarnings('ignore')

class SequenceCache:
    """
    Thread-safe sequence cache with preloaded BLOSUM62 matrix
    """
    _instance = None
    _lock = multiprocessing.Lock()
    
    def __new__(cls):
        with cls._lock:
            if cls._instance is None:
                cls._instance = super().__new__(cls)
                cls._instance.cache = {}
                cls._instance.aligner = PairwiseAligner()
                cls._instance.aligner.mode = 'global'
                cls._instance.aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
                cls._instance.aligner.open_gap_score = -11
                cls._instance.aligner.extend_gap_score = -1
            return cls._instance

    def get_sequences(self, filename):
        if filename not in self.cache:
            self.cache[filename] = self._read_sequences(filename)
        return self.cache[filename]
    
    def _read_sequences(self, filename):
        sequences = {}
        with open(filename, 'rt') as f:
            for record in SeqIO.parse(f, "fasta"):
                # Pre-compute sequence length and store with sequence
                seq_str = str(record.seq)
                sequences[record.id] = (seq_str, len(seq_str))
        return sequences

    def get_aligner(self):
        return self.aligner

@lru_cache(maxsize=10000)
def calculate_global_alignment_score(seq1, seq2):
    """
    Optimized alignment score calculation with early rejection
    """
    # Get cached aligner instance
    aligner = SequenceCache().get_aligner()
    
    # Quick rejection based on length difference
    len1, len2 = len(seq1), len(seq2)
    max_len = max(len1, len2)
    min_len = min(len1, len2)
    
    # More aggressive length-based filtering
    if min_len < max_len * 0.7:  # Sequences differ by more than 30%
        return 0.0
    
    # Quick composition check
    if abs(seq1.count('G') - seq2.count('G')) / max_len > 0.3:  # Example with glycine
        return 0.0
        
    actual_score = aligner.score(seq1, seq2)
    best_possible = max_len * 4
    worst_possible = -4 * max_len
    
    normalized_score = (actual_score - worst_possible) / (best_possible - worst_possible)
    return max(0.0, min(1.0, normalized_score))

def chunk_generator(total_pairs, chunk_size):
    """
    Generator for creating chunks to reduce memory usage
    """
    for start in range(0, total_pairs, chunk_size):
        end = min(start + chunk_size, total_pairs)
        yield start, end

def process_sequence_chunk(args):
    """
    Optimized chunk processing with early rejections
    """
    file1_path, file2_path, start_idx, end_idx = args
    
    results = []
    batch_size = 5000  # Increased batch size
    current_batch = []
    
    try:
        cache = SequenceCache()
        sequences1 = cache.get_sequences(file1_path)
        sequences2 = cache.get_sequences(file2_path)
        
        file1 = os.path.basename(file1_path)
        file2 = os.path.basename(file2_path)
        
        pairs_done = 0
        for seq1_id, (seq1, len1) in sequences1.items():
            for seq2_id, (seq2, len2) in sequences2.items():
                if start_idx <= pairs_done < end_idx:
                    score = calculate_global_alignment_score(seq1, seq2)
                    if score > 0.0:  # Only store non-zero scores
                        current_batch.append({
                            'File1': file1,
                            'File2': file2,
                            'Sequence1': seq1_id,
                            'Sequence2': seq2_id,
                            'Similarity_Score': score
                        })
                    
                    if len(current_batch) >= batch_size:
                        results.extend(current_batch)
                        current_batch = []
                
                pairs_done += 1
                if pairs_done >= end_idx:
                    break
            if pairs_done >= end_idx:
                break
        
        if current_batch:
            results.extend(current_batch)
        
        return results
    
    except Exception as e:
        print(f"Error processing chunk: {str(e)}")
        return []

def parallel_write_results(queue, output_file, chunk_size=10000):
    """
    Optimized batch writing of results with compression
    """
    buffer = []
    while True:
        results = queue.get()
        if results is None:
            if buffer:  # Write remaining results
                df = pd.DataFrame(buffer)
                df.to_csv(output_file, mode='a', header=False, index=False, compression='gzip')
            break
            
        buffer.extend(results)
        if len(buffer) >= chunk_size:
            df = pd.DataFrame(buffer)
            df.to_csv(output_file, mode='a', header=False, index=False, compression='gzip')
            buffer = []

def perform_all_vs_all_alignment(input_dir, output_dir, max_workers=None, chunk_size=50000):
    """
    Optimized alignment with larger chunks and better process management
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if max_workers is None:
        max_workers = max(1, multiprocessing.cpu_count() - 1)
    
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    detailed_output = os.path.join(output_dir, f'detailed_alignments_{timestamp}.csv.gz')
    summary_output = os.path.join(output_dir, f'similarity_matrix_{timestamp}.csv.gz')
    
    # Initialize output file with compression
    pd.DataFrame(columns=['File1', 'File2', 'Sequence1', 'Sequence2', 'Similarity_Score'])\
        .to_csv(detailed_output, index=False, compression='gzip')
    
    result_queue = Queue(maxsize=1000)  # Increased queue size
    writer_thread = Thread(target=parallel_write_results, 
                         args=(result_queue, detailed_output, 10000))  # Increased write batch size
    writer_thread.start()
    
    fasta_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.fasta')]
    file_pairs = list(combinations(fasta_files, 2))
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for pair_idx, (file1_path, file2_path) in enumerate(file_pairs, 1):
            try:
                # Pre-load sequences into cache
                cache = SequenceCache()
                cache.get_sequences(file1_path)
                cache.get_sequences(file2_path)
                
                seq_count1 = len(cache.get_sequences(file1_path))
                seq_count2 = len(cache.get_sequences(file2_path))
                total_pairs = seq_count1 * seq_count2
                
                print(f"\nProcessing pair {pair_idx}/{len(file_pairs)}: "
                      f"{os.path.basename(file1_path)} vs {os.path.basename(file2_path)}")
                
                # Create chunks using generator with larger chunk size
                chunks = [(file1_path, file2_path, start, end) 
                         for start, end in chunk_generator(total_pairs, chunk_size)]
                
                # Submit all chunks for this pair
                pair_futures = [executor.submit(process_sequence_chunk, chunk) for chunk in chunks]
                futures.extend(pair_futures)
                
            except Exception as e:
                print(f"Error processing file pair: {str(e)}")
                continue
        
        # Process all futures with progress bar
        for future in tqdm(futures, desc="Processing chunks"):
            results = future.result()
            if results:
                result_queue.put(results)
    
    result_queue.put(None)
    writer_thread.join()
    
    print("\nGenerating summary matrix...")
    try:
        summary_data = []
        for chunk in pd.read_csv(detailed_output, compression='gzip', chunksize=100000):
            summary = chunk.groupby(['File1', 'File2'])['Similarity_Score'].mean()
            summary_data.append(summary)
        
        final_summary = pd.concat(summary_data).groupby(level=[0, 1]).mean()
        final_summary.to_csv(summary_output, compression='gzip')
    except Exception as e:
        print(f"Error generating summary matrix: {str(e)}")
    
    return detailed_output, summary_output

def validate_directories(input_dir, output_dir):
    """
    Validate input and output directories
    """
    if not os.path.exists(input_dir):
        raise ValueError(f"Input directory does not exist: {input_dir}")
    
    fasta_files = [f for f in os.listdir(input_dir) if f.endswith('.fasta')]
    if not fasta_files:
        raise ValueError(f"No FASTA files found in input directory: {input_dir}")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    return True

def main():
    """
    Main function to run the sequence alignment pipeline
    """
    parser = argparse.ArgumentParser(
        description='Perform all-vs-all protein sequence similarity analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-i', '--input_dir', required=True,
                      help='Directory containing FASTA files')
    parser.add_argument('-o', '--output_dir', required=True,
                      help='Directory for output files')
    parser.add_argument('-w', '--workers', type=int, default=None,
                      help='Number of worker processes (default: CPU cores - 1)')
    parser.add_argument('-c', '--chunk_size', type=int, default=50000,
                      help='Size of chunks for parallel processing')
    parser.add_argument('--verbose', action='store_true',
                      help='Enable verbose output')
    
    args = parser.parse_args()
    
    try:
        # Validate directories
        validate_directories(args.input_dir, args.output_dir)
        
        # Set up logging level
        if args.verbose:
            print("Running in verbose mode")
        
        print(f"Starting analysis with {args.workers if args.workers else multiprocessing.cpu_count() - 1} workers")
        print(f"Input directory: {args.input_dir}")
        print(f"Output directory: {args.output_dir}")
        
        # Run the alignment
        start_time = datetime.datetime.now()
        detailed_file, summary_file = perform_all_vs_all_alignment(
            args.input_dir,
            args.output_dir,
            args.workers,
            args.chunk_size
        )
        end_time = datetime.datetime.now()
        
        # Print results
        print(f"\nAnalysis complete!")
        print(f"Time taken: {end_time - start_time}")
        print(f"Detailed results saved to: {detailed_file}")
        print(f"Summary matrix saved to: {summary_file}")
        
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    multiprocessing.freeze_support()  # For Windows compatibility
    main()