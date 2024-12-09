import sys
import pandas as pd
import os
import math
from multiprocessing import Pool, cpu_count

# Function to parse BLAST results
def parse_blast_results(file_path, chunksize=10000):
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    # Read the file with space delimiter in chunks
    return pd.read_csv(file_path, header=None, names=columns, skiprows=1, engine='python', chunksize=chunksize)

def calculate_composite_score(args):
    row, max_bit_score = args
    print(row)
    e_value = row['evalue']
    bit_score = row['bitscore']
    alignment_length = row['length']
    query_length = row['qend'] - row['qstart'] + 1
    percent_identity = row['pident']

    if e_value == 0:
        e_value = 1e-180  # Replace zero to avoid math errors

    # Normalized bit score
    normalized_bit_score = bit_score / max_bit_score
    log_e_value = -math.log10(e_value) if e_value > 0 else 0
    score_per_length = bit_score / alignment_length
    score_per_query_length = bit_score / query_length
    percent_identity_weight = percent_identity / 100
    coverage = alignment_length / query_length

    composite_score = (
        normalized_bit_score * log_e_value * score_per_length * 
        score_per_query_length * percent_identity_weight * coverage
    )

    return composite_score

def process_chunk(chunk, max_bit_score):
    with Pool(processes=cpu_count()) as pool:
        composite_scores = pool.map(calculate_composite_score, [(row, max_bit_score) for _, row in chunk.iterrows()])
    chunk['composite_score'] = composite_scores
    return chunk

def normalize_scores(file_path, chunksize=10000):
    chunks = parse_blast_results(file_path, chunksize)
    normalized_chunks = []
    
    for chunk in chunks:
        max_bit_score = chunk['bitscore'].max()
        chunk = process_chunk(chunk, max_bit_score)
        
        min_score = chunk['composite_score'].min()
        max_score = chunk['composite_score'].max()

        # Normalize the composite score to range between 0 and 1
        chunk['normalized_composite_score'] = (
            (chunk['composite_score'] - min_score) / (max_score - min_score)
        ).fillna(0)  # Handle division by zero if all scores are the same
        
        normalized_chunks.append(chunk)
    
    return pd.concat(normalized_chunks)

# Main script
if len(sys.argv) != 2:
    print("Usage: python normalize_blast_results.py <blast_csv_file>")
    sys.exit(1)

blast_output_file = sys.argv[1]  # Get file path from the command line argument
normalized_blast_df = normalize_scores(blast_output_file)

# Save the normalized results to a CSV file in the same directory as the input file
output_file = os.path.join(os.path.dirname(blast_output_file), 'normalized_blast_scores.csv')
normalized_blast_df.to_csv(output_file, sep='\t', index=False)
print(f"Normalized BLAST results saved to {output_file}")