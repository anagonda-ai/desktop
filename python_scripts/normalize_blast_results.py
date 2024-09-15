import pandas as pd
import os
import math

# Function to parse BLAST results
def parse_blast_results(file_path):
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    # Read the file with space delimiter
    blast_df = pd.read_csv(file_path, sep='\s+', header=None, names=columns, engine='python')
    return blast_df

def calculate_composite_score(row, max_bit_score):
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

def normalize_scores(df):
    # Calculate max and min composite scores
    df['composite_score'] = df.apply(lambda row: calculate_composite_score(row, df['bitscore'].max()), axis=1)
    min_score = df['composite_score'].min()
    max_score = df['composite_score'].max()

    # Normalize the composite score to range between 0 and 1
    df['normalized_composite_score'] = (
        (df['composite_score'] - min_score) / (max_score - min_score)
    ).fillna(0)  # Handle division by zero if all scores are the same

    return df


# Main script
blast_output_file = input("Enter path to BLAST CSV file: ")
blast_df = parse_blast_results(blast_output_file)
normalized_blast_df = normalize_scores(blast_df)

# Save the normalized results to a CSV file in the same directory as the input file
output_file = os.path.join(os.path.dirname(blast_output_file), 'normalized_blast_scores.csv')
normalized_blast_df.to_csv(output_file, sep='\t', index=False)
print(f"Normalized BLAST results saved to {output_file}")
