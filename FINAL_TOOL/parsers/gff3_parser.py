"""
GFF3 Parser
Parses GFF3 files and extracts gene information, ordering by position considering strand.
"""

import pandas as pd
from pathlib import Path
import logging
import re
from typing import List, Dict, Optional

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from FINAL_TOOL.utils.sequence_utils import load_genome_fasta, extract_gene_sequence

logger = logging.getLogger(__name__)


def parse_gff3_attributes(attributes_str: str) -> Dict[str, str]:
    """
    Parse GFF3 attributes column (9th column).
    
    Format: key1=value1;key2=value2;...
    
    Args:
        attributes_str: Attributes string from GFF3
        
    Returns:
        dict: Dictionary of attributes
    """
    attrs = {}
    if pd.isna(attributes_str) or not attributes_str:
        return attrs
    
    for item in str(attributes_str).split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key.strip()] = value.strip()
    
    return attrs


def extract_gene_id(attributes: Dict[str, str]) -> Optional[str]:
    """
    Extract gene ID from attributes.
    Tries multiple common attribute keys.
    """
    for key in ['ID', 'gene_id', 'GeneID', 'gene', 'Name']:
        if key in attributes:
            # Remove any prefix like "gene:" or "Gene:"
            gene_id = attributes[key]
            gene_id = re.sub(r'^gene:', '', gene_id, flags=re.IGNORECASE)
            gene_id = re.sub(r'^Gene:', '', gene_id, flags=re.IGNORECASE)
            return gene_id
    return None


def parse_gff3(gff3_file: str) -> pd.DataFrame:
    """
    Parse GFF3 file and extract gene features.
    
    Args:
        gff3_file: Path to GFF3 file
        
    Returns:
        DataFrame with columns: id, chromosome, start, end, strand, attributes
    """
    logger.info(f"Parsing GFF3 file: {gff3_file}")
    
    # GFF3 format columns
    columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    
    # Read GFF3 file, skipping comment lines
    data = []
    with open(gff3_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            
            # Only process gene features
            if parts[2].lower() == 'gene':
                data.append(parts[:9])
    
    if not data:
        logger.warning(f"No gene features found in {gff3_file}")
        return pd.DataFrame(columns=['id', 'chromosome', 'start', 'end', 'strand'])
    
    df = pd.DataFrame(data, columns=columns)
    
    # Convert coordinates to integers
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    
    # Extract gene IDs from attributes
    gene_ids = []
    for attr_str in df['attributes']:
        attrs = parse_gff3_attributes(attr_str)
        gene_id = extract_gene_id(attrs)
        if gene_id:
            gene_ids.append(gene_id)
        else:
            logger.warning(f"Could not extract gene ID from attributes: {attr_str[:100]}")
            gene_ids.append(None)
    
    df['id'] = gene_ids
    
    # Filter out rows without gene IDs
    df = df[df['id'].notna()].copy()
    
    # Create output dataframe
    result_df = pd.DataFrame({
        'id': df['id'],
        'chromosome': df['seqid'],
        'start': df['start'],
        'end': df['end'],
        'strand': df['strand']
    })
    
    logger.info(f"Parsed {len(result_df)} genes from GFF3 file")
    return result_df


def order_genes_by_position(df: pd.DataFrame) -> pd.DataFrame:
    """
    Order genes by chromosome and position, considering strand.
    
    For + strand: order by start (ascending)
    For - strand: order by start (descending) - genes on reverse strand
                   are ordered from end to start
    
    Args:
        df: DataFrame with columns: id, chromosome, start, end, strand
        
    Returns:
        Ordered DataFrame
    """
    df = df.copy()
    
    # Function to extract numeric part from chromosome for sorting
    def extract_numeric_part(chromosome):
        chromosome = str(chromosome)
        # Extract numeric part
        numeric = ''.join(filter(str.isdigit, chromosome))
        return int(numeric) if numeric else 0
    
    df['chromosome_numeric'] = df['chromosome'].apply(extract_numeric_part)
    
    # Sort by chromosome, then by strand, then by position
    # For + strand: ascending by start
    # For - strand: descending by start (to represent reverse order)
    def sort_key(row):
        chrom_num = row['chromosome_numeric']
        strand = str(row['strand'])
        # For reverse strand, use negative start to reverse order
        start = row['start'] if strand == '+' or strand == '1' else -row['start']
        return (chrom_num, row['chromosome'], start)
    
    df = df.sort_values(by=['chromosome_numeric', 'chromosome', 'start'])
    
    # For reverse strand genes, we want them in reverse order on their chromosome
    # Group by chromosome and strand, then reverse order for - strand
    ordered_rows = []
    for (chrom, strand), group in df.groupby(['chromosome', 'strand']):
        if strand == '-' or strand == '-1':
            # Reverse order for negative strand
            group = group.sort_values('start', ascending=False)
        else:
            # Normal order for positive strand
            group = group.sort_values('start', ascending=True)
        ordered_rows.append(group)
    
    if ordered_rows:
        result_df = pd.concat(ordered_rows, ignore_index=True)
    else:
        result_df = df
    
    result_df = result_df.drop(columns=['chromosome_numeric'])
    
    logger.info(f"Ordered {len(result_df)} genes by position")
    return result_df


def extract_sequences_from_fasta(gff3_file: str, fasta_file: str, output_csv: Optional[str] = None) -> pd.DataFrame:
    """
    Parse GFF3, order genes, and extract sequences from FASTA.
    
    Args:
        gff3_file: Path to GFF3 file
        fasta_file: Path to genome FASTA file
        output_csv: Optional path to save results CSV
        
    Returns:
        DataFrame with columns: id, chromosome, start, end, strand, sequence
    """
    # Parse GFF3
    genes_df = parse_gff3(gff3_file)
    
    if genes_df.empty:
        logger.warning("No genes found in GFF3 file")
        return pd.DataFrame(columns=['id', 'chromosome', 'start', 'end', 'strand', 'sequence'])
    
    # Order genes
    genes_df = order_genes_by_position(genes_df)
    
    # Load genome
    genome = load_genome_fasta(fasta_file)
    
    if not genome:
        logger.error(f"Failed to load genome from {fasta_file}")
        return genes_df
    
    # Extract sequences
    sequences = []
    for _, row in genes_df.iterrows():
        seq = extract_gene_sequence(
            genome,
            row['chromosome'],
            row['start'],
            row['end'],
            row['strand']
        )
        sequences.append(seq)
    
    genes_df['sequence'] = sequences
    
    # Filter out genes without sequences
    genes_df = genes_df[genes_df['sequence'].notna()].copy()
    
    logger.info(f"Extracted sequences for {len(genes_df)} genes")
    
    # Save if requested
    if output_csv:
        genes_df.to_csv(output_csv, index=False)
        logger.info(f"Saved results to {output_csv}")
    
    return genes_df

