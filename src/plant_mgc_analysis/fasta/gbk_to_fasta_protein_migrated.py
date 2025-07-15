
class GbkToFastaProteinProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: gbk_to_fasta_protein.py."""
    
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

# Output file for combined FASTA protein sequences
output_file = "./MGCs/combined_gene_proteins.fasta"

# Open output file for writing
with open(output_file, "w") as output_fh:
    # Dictionary to store the count of cluster names
    cluster_count = {}
    
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
                        # Extract the nucleotide sequence
                        nucleotide_seq = feature.extract(record.seq)
                        # Translate to protein sequence
                        protein_seq = nucleotide_seq.translate()
                        # Get gene ID (locus tag)
                        gene_id = feature.qualifiers.get("locus_tag", ["Unknown"])[0]
                        
                        # Increment count for cluster name
                        if cluster_name in cluster_count:
                            cluster_count[cluster_name] += 1
                        else:
                            cluster_count[cluster_name] = 1
                        
                        header = f">{cluster_name}_{cluster_count[cluster_name]}|{gene_id}|{organism_name}\n"
                        
                        output_fh.write(header)
                        output_fh.write(str(protein_seq) + "\n")

logger.info("Combined FASTA protein file created:", output_file)
