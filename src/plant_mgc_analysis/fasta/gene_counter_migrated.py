
class GeneCounterProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: gene_counter.py."""
    
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

def create_id_file_mapping(base_dir):
  id_file_mapping = {}
  
  for dirpath, _, filenames in os.walk(base_dir):
    # Exclude directories containing the ignore path (optional)
      for filename in filenames:
        if filename == "sequences_for_each_port_id.txt":
          file_path = os.path.join(dirpath, filename)
          
          with open(file_path, "r") as f:
            for line in f:
              if line.startswith('>'):
                protein_id = line.strip()[1:]
                id_file_mapping.setdefault(protein_id, []).append(file_path)
  
  return id_file_mapping

def count_total_sequences(base_dir):
  total_sequences = 0
  
  for dirpath, _, filenames in os.walk(base_dir):
    # Exclude directories containing the ignore path (optional)
    # if "/groups/itay_mayrose/alongonda/datasets/plantcyc/datasets/plantcyc" not in dirpath:
      for filename in filenames:
        if filename == "sequences_for_each_port_id.txt":
          # Build the file path
          file_path = os.path.join(dirpath, filename)
          
          # Count lines in the file
          try:
            with open(file_path, "r") as f:
              total_sequences += (sum(1 for line in f)) / 2 # 2 rows for each ID. One is the ID, the other is the sequence
            logger.info(f"Counted sequences in: {file_path}")
          except FileNotFoundError:
            logger.info(f"File not found: {file_path}")
          except Exception as e:
            logger.info(f"Error processing file: {file_path}, {e}")
  
  # Print the total number of sequences
  logger.info(f"Total number of sequences across all files: {total_sequences}")

if __name__ == "__main__":
    # Replace with your base directory
    base_dir = "/groups/itay_mayrose/alongonda/datasets/plantcyc/"
    count_total_sequences(base_dir)
    id_file_mapping = create_id_file_mapping(base_dir)
    plantcyc_in_others = 0
    for key, item in id_file_mapping.items():
        if "/groups/itay_mayrose/alongonda/datasets/plantcyc/datasets/plantcyc/15.1.1/data/sequences_for_each_port_id.txt" in item:
            plantcyc_in_others += 1
        logger.info(f"{key}, {item}")
    logger.info(f"{plantcyc_in_others} genes are in plantcyc and other dbs")
