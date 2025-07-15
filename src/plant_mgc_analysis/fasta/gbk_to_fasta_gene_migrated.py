
class GbkToFastaGeneProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: gbk_to_fasta_gene.py."""
    
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


#!/usr/bin/env python

from ..config.settings import get_settings
from ..core.base import BioinformaticsProcessor
from ..core.exceptions import ProcessingError, ValidationError
from ..utils.file_operations import file_manager
from loguru import logger
from Bio import SeqIO
import sys
import os

# Input directory containing GenBank files
input_directory = "./MGCs/"

# Output file for combined FASTA sequences
output_file = "./MGCs/combined_gene_split.fasta"

# Open output file for writing
with open(output_file, "w") as output_fh:
    # Iterate over GenBank files in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith(".gbk"):
            filepath = os.path.join(input_directory, filename)
            # Parse GenBank file
            for record in SeqIO.parse(filepath, "genbank"):
                # Get cluster name and organism name
                cluster_name = record.annotations["accessions"][0] if "accessions" in record.annotations else "UnknownCluster"
                organism_name = record.annotations["organism"] if "organism" in record.annotations else "UnknownOrganism"
                # Iterate over features in the record
                for feature in record.features:
                    if feature.type == "CDS":
                        # Extract the gene sequence
                        gene_seq = feature.extract(record.seq)
                        # Get gene ID (locus tag)
                        gene_id = feature.qualifiers.get("locus_tag", ["Unknown"])[0]
                        # Write gene sequence in FASTA format with cluster and organism names
                        header = f">{cluster_name}|{organism_name}|{gene_id}\n"
                        output_fh.write(header)
                        output_fh.write(str(gene_seq) + "\n")

logger.info("Combined FASTA file created:", output_file)
