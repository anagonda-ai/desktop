
class FastaSplitterProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: fasta_splitter.py."""
    
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

def split_fasta_by_gene(fasta_path, out_dir):
    basename = os.path.basename(fasta_path)
    gene_index = 1
    cluster_index = 0
    if "MGC_CANDIDATE" not in basename:
        gene_index = 0
        cluster_index = 1
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    with open(fasta_path, 'r') as f:
        gene = None
        cluster = None
        file_index = 1
        seq_lines = []
        for line in f:
            if line.startswith('>'):
                if gene:
                    # Write previous gene
                    cluster_dir = os.path.join(out_dir, cluster)
                    if not os.path.exists(cluster_dir):
                        os.makedirs(cluster_dir)
                    gene_file = os.path.join(cluster_dir, f"gene_{file_index}.fasta")
                    with open(gene_file, 'w') as gf:
                        gf.write(header)
                        gf.writelines(seq_lines)
                    file_index += 1
                # Start new gene
                header = line
                gene = line[1:].split('|')[gene_index].strip().replace(' ', '_')  # Use second field as gene name
                cluster = line[1:].split('|')[cluster_index].strip().replace('>', '')
                seq_lines = []
            else:
                seq_lines.append(line)
        # Write last gene
        if gene:
            cluster_dir = os.path.join(out_dir, cluster)
            if not os.path.exists(cluster_dir):
                os.makedirs(cluster_dir)
            gene_file = os.path.join(cluster_dir, f"gene_{file_index}.fasta")
            with open(gene_file, 'w') as gf:
                gf.write(header)
                gf.writelines(seq_lines)
            file_index += 1

def main(txt_file):
    txt_file_dir = os.path.dirname(txt_file)
    out_dir = os.path.join(txt_file_dir, "mgc_candidates_dir_fixed")
    with open(txt_file, 'r') as f:
        fasta_files = [line.strip() for line in f if line.strip()]
    for fasta_path in fasta_files:
        split_fasta_by_gene(fasta_path, out_dir)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        logger.info(f"Usage: python {sys.argv[0]} <path_to_txt_file>")
        sys.exit(1)
    main(sys.argv[1])