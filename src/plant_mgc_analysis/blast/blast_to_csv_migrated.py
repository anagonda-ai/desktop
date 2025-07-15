
class BlastToCsvProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: blast_to_csv.py."""
    
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

def convert_blast_output_to_csv(blast_output_file, csv_output_file):
    # Define the column headers
    headers = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]
    
    # Read the BLAST output file into a DataFrame
    df = pd.read_csv(blast_output_file, sep='\t', header=None, names=headers)
    
    # Convert relevant columns to float, coercing errors to NaN
    float_columns = ['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                     'sstart', 'send', 'evalue', 'bitscore']
    df[float_columns] = df[float_columns].apply(pd.to_numeric, errors='coerce')
    
    
    # Save the DataFrame to a CSV file
    df.to_csv(csv_output_file,index=False)

def main():
    blast_output_file = sys.argv[1]
    csv_output_file = sys.argv[2]
    
    convert_blast_output_to_csv(blast_output_file, csv_output_file)
    logger.info(f"Converted BLAST output to CSV: {csv_output_file}")

if __name__ == '__main__':
    main()