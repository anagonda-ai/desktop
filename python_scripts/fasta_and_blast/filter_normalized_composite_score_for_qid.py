import pandas as pd
import os

def filter_best_scores(df):
    # Group by 'qseqid' and keep the row with the highest 'normalized_composite_score' for each group
    idx = df.groupby('qseqid')['normalized_composite_score'].idxmax()
    return df.loc[idx]

# Main script
normalized_file_path = '/groups/itay_mayrose/alongonda/desktop/arabidopsis/normalized_blast_scores.csv'

# Load the normalized BLAST results
normalized_blast_df = pd.read_csv(normalized_file_path, sep='\t')

# Filter to keep only the row with the highest normalized_composite_score for each qseqid
best_scores_df = filter_best_scores(normalized_blast_df)

# Save the filtered results to a CSV file in the same directory as the input file
output_file = os.path.join(os.path.dirname(normalized_file_path), 'best_normalized_blast_scores.csv')
best_scores_df.to_csv(output_file, sep='\t', index=False)
print(f"Best normalized BLAST results saved to {output_file}")