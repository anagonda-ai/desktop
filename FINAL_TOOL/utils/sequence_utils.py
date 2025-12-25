"""
Utility functions for sequence operations
"""

import os
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


def load_genome_fasta(fasta_file):
    """
    Load genome FASTA file into a dictionary.
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        dict: Dictionary mapping chromosome/contig ID to Seq object
    """
    genome = {}
    
    try:
        with open(fasta_file, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                chrom_id = record.id.split()[0]  # Use first part of header as chromosome ID
                genome[chrom_id] = record.seq
                # Also try storing without any suffix
                if '.' in chrom_id:
                    genome[chrom_id.split('.')[0]] = record.seq
        
        logger.info(f"Loaded {len(genome)} chromosomes/contigs from {fasta_file}")
        return genome
    except Exception as e:
        logger.error(f"Error loading genome FASTA {fasta_file}: {e}")
        return {}


def extract_gene_sequence(genome, chromosome, start, end, strand='+'):
    """
    Extract gene sequence from genome using coordinates.
    
    Args:
        genome: Dictionary mapping chromosome ID to Seq object
        chromosome: Chromosome/contig identifier
        start: Start position (1-based, inclusive)
        end: End position (1-based, inclusive)
        strand: Strand ('+' or '-')
        
    Returns:
        str: Gene sequence (forward strand)
    """
    # Try different chromosome ID formats
    chrom_variants = [
        chromosome,
        str(chromosome),
        chromosome.lstrip('0'),
        f"chr{chromosome}",
        chromosome.split('.')[0]
    ]
    
    seq = None
    used_chrom = None
    
    for chrom_var in chrom_variants:
        if chrom_var in genome:
            seq = genome[chrom_var]
            used_chrom = chrom_var
            break
    
    if seq is None:
        logger.warning(f"Chromosome {chromosome} not found in genome")
        return None
    
    # Convert to 0-based indexing
    start_idx = start - 1
    end_idx = end
    
    # Extract sequence
    if start_idx < 0:
        start_idx = 0
    if end_idx > len(seq):
        end_idx = len(seq)
    
    gene_seq = seq[start_idx:end_idx]
    
    # Handle strand
    if strand == '-' or strand == '-1':
        gene_seq = gene_seq.reverse_complement()
    
    return str(gene_seq)


def extract_promoter_sequence(genome, chromosome, start, end, strand='+', upstream=1000):
    """
    Extract promoter sequence (upstream region) from genome.
    
    Args:
        genome: Dictionary mapping chromosome ID to Seq object
        chromosome: Chromosome/contig identifier
        start: Gene start position (1-based)
        end: Gene end position (1-based)
        strand: Strand ('+' or '-')
        upstream: Number of bp upstream to extract
        
    Returns:
        str: Promoter sequence
    """
    chrom_variants = [
        chromosome,
        str(chromosome),
        chromosome.lstrip('0'),
        f"chr{chromosome}",
        chromosome.split('.')[0]
    ]
    
    seq = None
    for chrom_var in chrom_variants:
        if chrom_var in genome:
            seq = genome[chrom_var]
            break
    
    if seq is None:
        logger.warning(f"Chromosome {chromosome} not found in genome")
        return None
    
    if strand == '+' or strand == '1':
        # Upstream is before start
        promoter_start = max(0, start - 1 - upstream)
        promoter_end = start - 1
        promoter_seq = seq[promoter_start:promoter_end]
    else:
        # Upstream is after end (reverse complement)
        promoter_start = end
        promoter_end = min(len(seq), end + upstream)
        promoter_seq = seq[promoter_start:promoter_end].reverse_complement()
    
    return str(promoter_seq)


