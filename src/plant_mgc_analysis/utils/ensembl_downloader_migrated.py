
class EnsemblDownloaderProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: ensembl_downloader.py."""
    
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
import ftplib
import os

def download_cdna_folders(ftp, path, download_dir):
    try:
        items = ftp.nlst(path)
        for item in items:
            species_name = item.split('/')[-1]
            cdna_path = os.path.join(item, 'pep/')
            try:
                ftp.cwd(cdna_path)
                logger.info(f"Found cdna directory for {species_name}, downloading...")
                species_download_path = os.path.join(download_dir, species_name, 'pep')
                os.makedirs(species_download_path, exist_ok=True)
                cdna_files = ftp.nlst(cdna_path)
                for cdna_file in cdna_files:
                    local_file_path = os.path.join(species_download_path, os.path.basename(cdna_file))
                    with open(local_file_path, 'wb') as local_file:
                        ftp.retrbinary(f'RETR {cdna_file}', local_file.write)
                        logger.info(f"Downloaded {local_file_path}")
            except ftplib.error_perm:
                logger.info(f"No cdna directory found for {species_name}.")
            finally:
                ftp.cwd('..')
    except Exception as e:
        logger.info(f"Error accessing {path}: {e}")

def main():
    # Get user input for the download directory
    download_dir = input("Enter the directory to save downloaded cdna files: ")
    os.makedirs(download_dir, exist_ok=True)

    ftp_url = 'ftp.ensemblgenomes.ebi.ac.uk'
    base_path = '/pub/plants/release-59/fasta/'

    ftp = ftplib.FTP(ftp_url)
    ftp.login()
    download_cdna_folders(ftp, base_path, download_dir)
    ftp.quit()

if __name__ == "__main__":
    main()
