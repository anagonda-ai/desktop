
class HeaderExtractorProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: header_extractor.py."""
    
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
import concurrent.futures
import csv
import os

def extract_info_from_header(header):
    # Regular expression to extract gene and transcript name, start, and end positions
    parts = header.split(' ')
    gene_name = parts[3].split(':')[1]
    transcript_name = parts[4].split(':')[1]
    chrosomose = parts[2].split(':')[2]
    start_position = int(parts[2].split(':')[3])
    end_position = int(parts[2].split(':')[4])
    strand = int(parts[2].split(':')[5])
    return gene_name, transcript_name, chrosomose, start_position, end_position, strand

def extract_numeric_part(chromosome):
    return int(''.join(filter(str.isdigit, chromosome)) or 0)

def process_file(file_path):
    data = {}
    with open(file_path, 'r') as f:
        current_header = None
        current_sequence = []

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_sequence:
                    gene_name, transcript_name, chrosomose, start_position, end_position, strand = extract_info_from_header(current_header)
                    protein_sequence = ''.join(current_sequence)
                    transcript_length = end_position - start_position

                    if gene_name not in data or transcript_length > data[gene_name][4] - data[gene_name][3]:
                        data[gene_name] = (gene_name, transcript_name, chrosomose, start_position, end_position, strand, protein_sequence)
                
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)
        
        if current_header and current_sequence:
            gene_name, transcript_name, chrosomose, start_position, end_position, strand = extract_info_from_header(current_header)
            protein_sequence = ''.join(current_sequence)
            transcript_length = end_position - start_position

            if gene_name not in data or transcript_length > data[gene_name][4] - data[gene_name][3]:
                data[gene_name] = (gene_name, transcript_name, chrosomose, start_position, end_position, strand, protein_sequence)
    
    sorted_data = sorted(data.values(), key=lambda x: (extract_numeric_part(x[2]), x[2], x[3], x[4]))

    output_file = os.path.join(os.path.dirname(file_path), 'extracted_data.csv')
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['gene_name', 'transcript_name', 'chrosomose', 'start_position', 'end_position', 'strand','protein_sequence'])
        writer.writerows(sorted_data)
    
    logger.info(f"Data saved to {output_file}")

def process_fa_files(root_dir):
    file_paths = []
    for subdir, _, files in os.walk(root_dir):
        for file in files:
            if file.endswith('.fasta') or file.endswith('.fa'):
                file_path = os.path.join(subdir, file)
                file_paths.append(file_path)
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = {executor.submit(process_file, file_path): file_path for file_path in file_paths}
        for future in concurrent.futures.as_completed(futures):
            file_path = futures[future]
            try:
                future.result()
                logger.info(f"Processed file: {file_path}")
            except Exception as exc:
                logger.info(f"Error processing file {file_path}: {exc}")

def main():
    root_dir = "/groups/itay_mayrose/alongonda/datasets/full_genomes/ensembl/organisms"
    process_fa_files(root_dir)

if __name__ == "__main__":
    main()