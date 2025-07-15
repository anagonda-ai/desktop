
class SubmitMultipleJobsProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: submit_multiple_jobs.py."""
    
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
import subprocess

import pandas as pd

full_genome_dir = "/groups/itay_mayrose/alongonda/datasets/full_genomes"
# Define genome directories

genome_dirs = [
    os.path.join(full_genome_dir, "ensembl/processed_annotations_test_no_chloroplast_with_sequences"),
    os.path.join(full_genome_dir, "phytozome/processed_annotations_with_chromosomes_no_chloroplast_with_sequences"),
    os.path.join(full_genome_dir, "plaza/processed_annotations_with_chromosomes_no_chloroplast_with_sequences")
]
kegg_db = "/groups/itay_mayrose/alongonda/datasets/KEGG/fasta/merged_metabolic_pathways/merged_metabolic_pathways_db"
head_output_dir = "/groups/itay_mayrose/alongonda/Plant_MGC/kegg_output"
output_dir = os.path.join(head_output_dir, "kegg_scanner_min_genes_based_metabolic")
temp_dir = os.path.join(head_output_dir, "blast_temp_annotated_metabolic")
annotated_dir = os.path.join(head_output_dir, "annotated_genomes_metabolic")

genome_files = []
for genome_dir in genome_dirs:
    for file in os.listdir(genome_dir):
        if file.endswith('.csv'):
            genome_files.append(os.path.join(genome_dir, file))
            
for genome_file in genome_files:
    blast_output = os.path.join(temp_dir, f"blast_result_{os.path.basename(genome_file)}.xml")
    if not os.path.exists(blast_output):
        
        df = pd.read_csv(genome_file)
        fasta_query = os.path.join(temp_dir, os.path.basename(genome_file).replace('.csv', '.fasta'))
        # Save genome sequences to temporary FASTA
        if not os.path.exists(fasta_query):
            with open(fasta_query, "w") as f:
                for idx, row in df.iterrows():
                    f.write(f">{row['id']}\n{row['sequence']}\n")
        filename = os.path.basename(genome_file).replace('.csv', '_annotated.csv')
        output_path = os.path.join(annotated_dir, filename)
        subprocess.run("sbatch /groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/blast_with_params.sh {} {} {}".format(
            fasta_query, kegg_db, blast_output), shell=True)