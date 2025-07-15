
class PlazaDownloaderProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: plaza_downloader.py."""
    
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
import requests
from bs4 import BeautifulSoup
import os

def download_fasta_files(parent_url, download_dir):
    # Create the download directory if it doesn't exist
    os.makedirs(download_dir, exist_ok=True)

    # Fetch the HTML content of the parent URL
    try:
        response = requests.get(parent_url)
        response.raise_for_status()  # Raise an error for bad responses
    except requests.exceptions.RequestException as e:
        logger.info(f"Error fetching the parent URL: {e}")
        return

    # Parse the HTML using BeautifulSoup
    soup = BeautifulSoup(response.content, 'html.parser')

    # Find the div with the ID 'fasta_data_genome'
    fasta_div = soup.find('div', id='fasta_data_genome')
    
    # Find all links ending with .fasta.gz within that div
    fasta_links = [a['href'] for a in fasta_div.find_all('a', href=True) if a['href'].endswith('.fasta.gz')]
    
    # Debug output to show found links
    logger.info(f"Found {len(fasta_links)} fasta.gz links.")

    if not fasta_links:
        logger.info("No fasta.gz files found.")
        return

    # Download each fasta.gz file
    for link in fasta_links:
        # Construct the full URL
        if not link.startswith('http'):
            link = f"{parent_url.rsplit('/', 1)[0]}/{link}"  # Ensure it's a full URL

        logger.info(f"Downloading {link}...")
        try:
            fasta_response = requests.get(link)
            fasta_response.raise_for_status()

            # Write the content to a file
            file_name = os.path.join(download_dir, os.path.basename(link))
            with open(file_name, 'wb') as file:
                file.write(fasta_response.content)
                logger.info(f"Downloaded {file_name}")
        except requests.exceptions.RequestException as e:
            logger.info(f"Error downloading {link}: {e}")

if __name__ == "__main__":
    # URL of the parent directory
    parent_url = "https://bioinformatics.psb.ugent.be/plaza/versions/plaza_v5_dicots/download/download"
    # Directory to save the downloaded files
    download_dir = input("Enter the directory to save downloaded fasta.gz files: ")

    download_fasta_files(parent_url, download_dir)
