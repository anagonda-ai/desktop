
class GbkToFastaProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: gbk_to_fasta.py."""
    
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


#!/usr/bin/env python

from ..config.settings import get_settings
from ..core.base import BioinformaticsProcessor
from ..core.exceptions import ProcessingError, ValidationError
from ..utils.file_operations import file_manager
from loguru import logger
from Bio import SeqIO
import sys
import os

# Input directory containing GenBank files
input_directory = "./MGCs/"

# Output file for combined FASTA sequences
output_file = "./MGCs/combined.fasta"

# Open output file for writing
with open(output_file, "w") as output_fh:
    # Iterate over GenBank files in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith(".gbk"):
            filepath = os.path.join(input_directory, filename)
            # Parse GenBank file and write sequences in FASTA format to output file
            for record in SeqIO.parse(filepath, "genbank"):
                SeqIO.write(record, output_fh, "fasta")

logger.info("Combined FASTA file created:", output_file)
