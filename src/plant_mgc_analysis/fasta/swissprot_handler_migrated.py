
class SwissprotHandlerProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: swissprot_handler.py."""
    
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
from Bio import ExPASy

def fetch_protein_sequence(swiss_prot_id):
    try:
        handle = ExPASy.get_sprot_raw(swiss_prot_id)
        data = handle.read()
        handle.close()

        # Find the line containing the sequence
        sequence_start_index = data.find("SQ   SEQUENCE   ") + len("SQ   SEQUENCE   ")
        sequence_end_index = data.find("//", sequence_start_index)
        sequence_lines = data[sequence_start_index:sequence_end_index].split("\n")

        # Extract the sequence lines and concatenate them
        protein_sequence = "".join(line.strip() for line in sequence_lines[1:])
        return protein_sequence
    except Exception as e:
        logger.info(f"Failed to retrieve protein sequence for ID: {swiss_prot_id}, error: {e}")
        return None

def process_swissprotids_file(swissprotids_file):
    output_file = os.path.join(os.path.dirname(swissprotids_file), "sequences_for_each_port_id.txt")
    
    with open(swissprotids_file, "r") as f:
        swiss_prot_ids = f.read().splitlines()

    with open(output_file, "w") as f:
        for swiss_prot_id in swiss_prot_ids:
            protein_sequence = fetch_protein_sequence(swiss_prot_id)
            if protein_sequence:
                f.write(f">{swiss_prot_id}\n")
                f.write(protein_sequence + "\n")
                logger.info(f"Retrieved protein sequence for ID: {swiss_prot_id}")
            else:
                logger.info(f"Failed to retrieve protein sequence for ID: {swiss_prot_id}")
    
    logger.info(f"Sequences saved in: {output_file}")

def find_swissprotids_files(base_dir):
    ignore_path = "/groups/itay_mayrose/alongonda/datasets/plantcyc/datasets/plantcyc"
    
    for dirpath, _, filenames in os.walk(base_dir):
        if ignore_path not in dirpath:  # Exclude directories containing the ignore path
            logger.info(dirpath)
            for filename in filenames:
                if filename == "swissprotids.txt":
                    swissprotids_file = os.path.join(dirpath, filename)
                    logger.info(f"Processing file: {swissprotids_file}")
                    process_swissprotids_file(swissprotids_file)

if __name__ == "__main__":
    base_dir = "/groups/itay_mayrose/alongonda/datasets/plantcyc/"  # Replace with your base directory
    find_swissprotids_files(base_dir)
