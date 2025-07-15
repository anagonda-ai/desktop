
class ColumnRenamerProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: column_renamer.py."""
    
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
import pandas as pd

def load_and_rename_csv(input_csv, output_csv, new_column_name):
    # Load the CSV file into a DataFrame
    df = pd.read_csv(input_csv)
    
    # Rename the first column
    first_column_name = df.columns[0]
    df.rename(columns={first_column_name: new_column_name}, inplace=True)
    
    # Save the modified DataFrame to a new CSV file
    df.to_csv(output_csv, index=False)
    logger.info(f"Renamed first column to '{new_column_name}' and saved to {output_csv}")

def main():
    input_csv = "/groups/itay_mayrose/alongonda/datasets/full_genomes/annotations/phytozome_annotations.csv"
    output_csv = "/groups/itay_mayrose/alongonda/datasets/full_genomes/annotations/phytozome_annotations_fixed.csv"
    new_column_name = "id"
    
    load_and_rename_csv(input_csv, output_csv, new_column_name)

if __name__ == "__main__":
    main()