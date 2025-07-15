
class FastaCreatorE2p2Processor(BioinformaticsProcessor):
    """Migrated from legacy script: fasta_creator_e2p2.py."""
    
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
import csv
import os
from Bio import SeqIO

# Input files
CSV_FILE = "/groups/itay_mayrose/alongonda/datasets/plantcyc/pmn_mgc_potential/mgc_candidates_process/negative_training_set.csv"
FASTA_FILE = "/groups/itay_mayrose/alongonda/datasets/plantcyc/pmn_mgc_potential/mgc_candidates_process/pmn.fasta"
OUTPUT_DIR = "/groups/itay_mayrose/alongonda/datasets/plantcyc/pmn_mgc_potential/mgc_candidates_process/negative_candidates_fasta_files"

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

def parse_gene_ids(gene_id_field):
    """Parse gene IDs from a comma-separated string."""
    return gene_id_field.split(",")

def read_csv(file_path):
    """Read the CSV file and parse the pathway-gene map."""
    entries = []
    negative_candidate_id = 1
    with open(file_path, "r") as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            pathway = row["Pathway (Occurrence)"].split(" (")[0]  # Remove occurrence info
            gene_ids = parse_gene_ids(row["Window Genes"])
            for gene_id in gene_ids:
                entries.append({
                    "negative_candidate_id": f"NEGATIVE_CANDIDATE_{negative_candidate_id}",
                    "gene_id": gene_id,
                    "pathway": pathway
                })
            negative_candidate_id += 1
    return entries

def load_fasta_records(fasta_file):
    """Load all FASTA records into a dictionary for fast lookup."""
    return {record.id.split(" ")[0]: record for record in SeqIO.parse(fasta_file, "fasta")}

def write_fasta_files(entries, fasta_records, output_dir):
    """Write individual FASTA files for each candidate."""
    for entry in entries:
        gene_id = entry["gene_id"]
        negative_candidate_id = entry["negative_candidate_id"]
        pathway = entry["pathway"]
        
        if gene_id in fasta_records:
            record = fasta_records[gene_id]
            organism = record.description.split(" | ")[1] if " | " in record.description else "unknown"
            new_header = f">{negative_candidate_id} | {gene_id} | {pathway} | {organism}"
            fasta_filename = os.path.join(output_dir, f"{negative_candidate_id}.fasta")
            
            with open(fasta_filename, "a") as fasta_file:
                fasta_file.write(f"{new_header}\n{str(record.seq)}\n")

def main():
    entries = read_csv(CSV_FILE)
    fasta_records = load_fasta_records(FASTA_FILE)
    write_fasta_files(entries, fasta_records, OUTPUT_DIR)
    logger.info(f"FASTA files created in folder: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()