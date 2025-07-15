
class FastaOrdererProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: fasta_orderer.py."""
    
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
import os
import sys
from Bio import SeqIO
import re

def parse_fasta_and_positions(fasta_file):
    gene_positions = []
    sequences = {}
    headers = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        # Extract information from the description
        description = record.description
        match = re.search(r'chromosome:TAIR10:(\d+):(\d+):(\d+):(-?\d+)', description)
        if match:
            chrom, start, end, strand = match.groups()
            position = (int(chrom), int(start), int(end), int(strand))
            gene_positions.append((record.id, position))
            sequences[record.id] = str(record.seq)
            headers[record.id] = description

    return gene_positions, sequences, headers

def sort_sequences_by_position(gene_positions, sequences):
    # Sort genes by their start position
    sorted_genes = sorted(gene_positions, key=lambda x: (x[1][0], x[1][1], x[1][2], x[1][3]))
    sorted_sequences = {gene_id: sequences[gene_id] for gene_id, _ in sorted_genes}
    return sorted_sequences

def write_ordered_fasta(output_file, sorted_sequences, headers):
    with open(output_file, 'w') as out_file:
        for gene_id, sequence in sorted_sequences.items():
            out_file.write(f">{headers[gene_id]}\n")
            out_file.write(f"{sequence}\n")    
            
def main():
    if len(sys.argv) != 2:
        logger.info("Usage: python order_fasta_by_location.py <input_fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_dir = os.path.dirname(fasta_file)
    output_file = os.path.join(output_dir, "sorted_" + os.path.basename(fasta_file))

    # Parse the file and get gene positions and sequences
    gene_positions, sequences, headers = parse_fasta_and_positions(fasta_file)

    # Sort sequences by their positions
    sorted_sequences = sort_sequences_by_position(gene_positions, sequences)

    # Write the sorted sequences to the output file
    write_ordered_fasta(output_file, sorted_sequences, headers)

if __name__ == "__main__":
    main()
