
class PhytozomeNameExtractorProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: phytozome_name_extractor.py."""
    
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
import re
import requests
import os
import pandas as pd

def extract_dataset_id(filename):
    """
    Extracts the dataset ID from a filename.
    The dataset ID is the numeric part before "_v" or right before ".gene_transformed_filtered".
    """
    match = re.search(r'_\d{3}_', filename)
    return match.group(0).strip('_') if match else None

def get_organism_name(dataset_id):
    """Retrieves the organism name from Phytozome using the correct API."""
    api_url = f"https://files.jgi.doe.gov/phytozome_file_list/?api_version=2&d=asc&q={dataset_id}"
    response = requests.get(api_url)

    if response.status_code == 200:
        try:
            data = response.json()
            if "organisms" in data and len(data["organisms"]) > 0:
                organism_info = data["organisms"][0]  # Extract first result
                try:
                    species = organism_info["top_hit"]["metadata"]["ncbi_taxon"]['ncbi_taxon_species']
                except KeyError:
                    species = organism_info["top_hit"]["metadata"]["phytozome"]["proteome_name"]
                return f"{species}"
        except ValueError:
            return "Invalid JSON response"
    return f"Failed request: {response.status_code}"

def process_datasets(directory):
    """
    Iterates through filenames in the dataset directory, extracts dataset IDs, and fetches organism names.
    """
    if not os.path.exists(directory):
        logger.info(f"Error: Directory '{directory}' does not exist.")
        return

    dataset_mapping = {}

    for filename in os.listdir(directory):
        logger.info(f"Processing: {filename}")
        file_path = os.path.join(directory, filename)
        if os.path.isdir(file_path):  # Ensure it's a subdir, not a file
            dataset_id = extract_dataset_id(filename)
            if dataset_id:
                organism_name = get_organism_name(dataset_id)
                dataset_mapping[filename] = (dataset_id, organism_name)
                logger.info(f"Processed: {filename} -> Dataset ID: {dataset_id}, Organism: {organism_name}")
    
    return dataset_mapping

if __name__ == "__main__":
    DATASETS_PATH = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/blast_results_chromosome_separated/best_hits_by_organism"
    dataset_info = process_datasets(DATASETS_PATH)

    # Save results to a file (optional)
    # Convert the dataset mapping to a DataFrame
    df = pd.DataFrame.from_dict(dataset_info, orient='index', columns=['Dataset ID', 'Organism'])
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'Original Filename'}, inplace=True)

    # Save the DataFrame to a CSV file
    outpud_dir = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/output"
    os.makedirs(outpud_dir, exist_ok=True)
    output_csv_file = "dataset_organism_mapping.csv"
    output_path = os.path.join(outpud_dir, output_csv_file)
    df.to_csv(output_path)

    logger.info(f"\nMapping saved to {output_path}")