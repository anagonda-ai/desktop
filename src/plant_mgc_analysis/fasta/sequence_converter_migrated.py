
class SequenceConverterProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: sequence_converter.py."""
    
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
from Bio.Seq import Seq
import os
from Bio import SeqIO

input_dir = "/groups/itay_mayrose_nosnap/ronenshtein/kalanchoe_2024/prep/run/pt_genes/"
output_dir = "/groups/itay_mayrose/alongonda/datasets/full_genomes/chloroplast_prot_genes/"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Iterate over all fasta files in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith(".fasta") or filename.endswith(".fa"):
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, f"translated_{filename}")

        with open(input_path, "r") as input_handle, open(output_path, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                # Translate the nucleotide sequence to a protein sequence
                # Create a Seq object
                dna_seq = Seq(record.seq)
                protein_seq = str(dna_seq.translate())
                # Write the translated sequence to the output file
                output_handle.write(f">{record.id}\n{protein_seq}\n")
                logger.info(f"Translated and saved {record.id} to {output_path}")
