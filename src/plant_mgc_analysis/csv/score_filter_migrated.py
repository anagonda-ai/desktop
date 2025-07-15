
class ScoreFilterProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: score_filter.py."""
    
    def __init__(self, **kwargs):
        """Initialize processor."""
        super().__init__(**kwargs)
        self.settings = get_settings()
    
    def validate_input(self, data):
        """Validate input data."""
        pass  # Implement validation
    
    def process(self, data, **kwargs):
        """Process data."""
        # Original script logic here
        pass


from ..config.settings import get_settings
from ..core.base import BioinformaticsProcessor
from ..core.exceptions import ProcessingError, ValidationError
from ..utils.file_operations import file_manager
from loguru import logger
import sys
import pandas as pd
import os

def filter_best_scores(df):
    # Group by 'qseqid' and keep the row with the highest 'normalized_composite_score' for each group
    idx = df.groupby('qseqid')['normalized_composite_score'].idxmax()
    return df.loc[idx]

# Main script
if len(sys.argv) != 2:
    logger.info("Usage: python filter_normalized_composite_score_for_qid.py <normalized_blast_csv_file>")
    sys.exit(1)

normalized_file_path = sys.argv[1]  # Get file path from the command line argument

# Load the normalized BLAST results
normalized_blast_df = pd.read_csv(normalized_file_path, sep='\t')

# Filter to keep only the row with the highest normalized_composite_score for each qseqid
best_scores_df = filter_best_scores(normalized_blast_df)

# Save the filtered results to a CSV file in the same directory as the input file
output_file = os.path.join(os.path.dirname(normalized_file_path), 'best_normalized_blast_scores.csv')
best_scores_df.to_csv(output_file, sep='\t', index=False)
logger.info(f"Best normalized BLAST results saved to {output_file}")