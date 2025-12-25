"""
KEGG Annotation Module
BLASTs gene sequences against KEGG metabolic database to identify metabolic genes.
"""

import os
import subprocess
import pandas as pd
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Tuple
import logging
import tempfile

from Bio.Blast import NCBIXML

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from FINAL_TOOL import config

logger = logging.getLogger(__name__)


def create_gene_fasta(genes_df: pd.DataFrame, output_fasta: str, id_mapping_file: str) -> None:
    """
    Create FASTA file from gene sequences with short IDs and mapping file.
    
    Args:
        genes_df: DataFrame with 'id' and 'sequence' columns
        output_fasta: Path to output FASTA file
        id_mapping_file: Path to output ID mapping file (short_id -> original_id)
    """
    logger.info(f"Creating gene FASTA file: {output_fasta}")
    
    with open(output_fasta, 'w') as fasta_file, open(id_mapping_file, 'w') as mapping_file:
        for i, row in genes_df.iterrows():
            gene_id = str(row['id'])
            seq = str(row.get('sequence', ''))
            
            if not seq or seq == 'nan':
                continue
            
            # Remove stop codons if present
            seq = seq.replace('*', '')
            
            # Create short ID for BLAST
            short_id = f"gene_{i+1:05d}"
            
            # Write FASTA
            fasta_file.write(f">{short_id}\n{seq}\n")
            
            # Write mapping
            mapping_file.write(f"{short_id}\t{gene_id}\n")
    
    logger.info(f"Created FASTA with {len(genes_df)} genes")


def create_blast_db(fasta_file: str, db_name: str = None) -> str:
    """
    Create BLAST database from FASTA file.
    
    Args:
        fasta_file: Path to FASTA file
        db_name: Optional database name (defaults to fasta_file without extension)
        
    Returns:
        Database name (without extension)
    """
    if db_name is None:
        db_name = str(Path(fasta_file).with_suffix(''))
    
    logger.info(f"Creating BLAST database: {db_name}")
    
    cmd = [
        'makeblastdb',
        '-in', fasta_file,
        '-dbtype', 'prot',
        '-parse_seqids',
        '-out', db_name
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info("BLAST database created successfully")
        return db_name
    except subprocess.CalledProcessError as e:
        logger.error(f"Error creating BLAST database: {e.stderr}")
        raise


def run_blastp(query_fasta: str, db_name: str, output_xml: str, 
               evalue: float = 1e-3, num_threads: int = 4, 
               qcov_hsp_perc: float = 70.0) -> None:
    """
    Run BLASTP from query FASTA against database.
    
    Args:
        query_fasta: Path to query FASTA file (KEGG)
        db_name: BLAST database name
        output_xml: Path to output XML file
        evalue: E-value threshold
        num_threads: Number of threads
        qcov_hsp_perc: Query coverage percentage
    """
    if os.path.exists(output_xml):
        logger.info(f"BLAST output already exists: {output_xml}")
        return
    
    logger.info(f"Running BLASTP: {query_fasta} vs {db_name}")
    
    cmd = [
        'blastp',
        '-task', 'blastp-fast',
        '-query', query_fasta,
        '-db', db_name,
        '-out', output_xml,
        '-outfmt', '5',  # XML format
        '-evalue', str(evalue),
        '-qcov_hsp_perc', str(qcov_hsp_perc),
        '-num_threads', str(num_threads),
        '-max_target_seqs', '5'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info("BLASTP completed successfully")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running BLASTP: {e.stderr}")
        raise


def parse_blast_xml(blast_xml: str, id_mapping_file: str, 
                    evalue_threshold: float = 1e-3,
                    identity_threshold: float = 50.0,
                    coverage_threshold: float = 70.0) -> Dict[str, List[Tuple[str, str, str]]]:
    """
    Parse BLAST XML results and extract KEGG annotations.
    
    Args:
        blast_xml: Path to BLAST XML output
        id_mapping_file: Path to ID mapping file (short_id -> original_id)
        evalue_threshold: E-value threshold
        identity_threshold: Identity percentage threshold
        coverage_threshold: Coverage percentage threshold
        
    Returns:
        Dictionary mapping gene_id to list of (kegg_id, annotation, pathway) tuples
    """
    logger.info(f"Parsing BLAST results: {blast_xml}")
    
    # Load ID mapping
    id_mapping = pd.read_csv(id_mapping_file, sep='\t', header=None, names=['short', 'orig'])
    id_dict = dict(zip(id_mapping['short'], id_mapping['orig']))
    
    # Dictionary to store gene -> KEGG mappings
    gene_to_kegg = defaultdict(list)
    
    with open(blast_xml, 'r') as handle:
        for record in NCBIXML.parse(handle):
            query_def = record.query
            query_len = record.query_length
            
            # KEGG entries have format: kegg_id$annotation$pathway
            if '$' not in query_def:
                continue
            
            parts = query_def.split('$')
            if len(parts) < 3:
                continue
            
            kegg_id, annotation, pathway = map(str.strip, parts[:3])
            
            # Find best hit for this KEGG entry
            best_hit = None
            best_score = -1
            
            for alignment in record.alignments:
                hit_id = alignment.hit_id.strip().split()[0]
                
                for hsp in alignment.hsps:
                    identity = 100 * hsp.identities / hsp.align_length
                    coverage = 100 * hsp.align_length / query_len
                    
                    if (hsp.expect < evalue_threshold and
                        identity >= identity_threshold and
                        coverage >= coverage_threshold):
                        if hsp.bits > best_score:
                            best_hit = hit_id
                            best_score = hsp.bits
                    break  # Only consider first HSP
            
            # Map back to original gene ID
            if best_hit and best_hit in id_dict:
                orig_id = id_dict[best_hit]
                gene_to_kegg[orig_id].append((kegg_id, annotation, pathway))
    
    logger.info(f"Found KEGG annotations for {len(gene_to_kegg)} genes")
    return dict(gene_to_kegg)


def annotate_genes_with_kegg(genes_df: pd.DataFrame, 
                             kegg_db: str = None,
                             temp_dir: str = None,
                             evalue: float = None,
                             identity: float = None,
                             coverage: float = None) -> pd.DataFrame:
    """
    Annotate genes with KEGG information using BLAST.
    
    Args:
        genes_df: DataFrame with 'id' and 'sequence' columns
        kegg_db: Path to KEGG FASTA database (default from config)
        temp_dir: Temporary directory for intermediate files (default from config)
        evalue: E-value threshold (default from config)
        identity: Identity threshold (default from config)
        coverage: Coverage threshold (default from config)
        
    Returns:
        DataFrame with added 'kegg_ids', 'annotation', and 'pathway' columns
    """
    # Use defaults from config
    if kegg_db is None:
        kegg_db = config.KEGG_DB
    if temp_dir is None:
        temp_dir = config.TEMP_DIR
    if evalue is None:
        evalue = config.BLAST_EVALUE
    if identity is None:
        identity = config.BLAST_IDENTITY
    if coverage is None:
        coverage = config.BLAST_COVERAGE
    
    # Create temp directory
    os.makedirs(temp_dir, exist_ok=True)
    
    # Create temporary file paths
    gene_fasta = os.path.join(temp_dir, "genes.fasta")
    id_mapping = os.path.join(temp_dir, "id_mapping.tsv")
    blast_db = os.path.join(temp_dir, "genes_db")
    blast_xml = os.path.join(temp_dir, "blast_results.xml")
    
    # Step 1: Create gene FASTA
    create_gene_fasta(genes_df, gene_fasta, id_mapping)
    
    # Step 2: Create BLAST database
    create_blast_db(gene_fasta, blast_db)
    
    # Step 3: Run BLASTP (KEGG as query, genes as database - reverse BLAST)
    run_blastp(kegg_db, blast_db, blast_xml, evalue=evalue, qcov_hsp_perc=coverage)
    
    # Step 4: Parse BLAST results
    gene_to_kegg = parse_blast_xml(blast_xml, id_mapping, 
                                   evalue_threshold=evalue,
                                   identity_threshold=identity,
                                   coverage_threshold=coverage)
    
    # Step 5: Add annotations to DataFrame
    result_df = genes_df.copy()
    
    kegg_ids_list = []
    annotations_list = []
    pathways_list = []
    
    for gene_id in result_df['id']:
        hits = gene_to_kegg.get(str(gene_id), [])
        if hits:
            kegg_ids_list.append(','.join(k[0] for k in hits))
            annotations_list.append(','.join(k[1] for k in hits))
            pathways_list.append(','.join(k[2] for k in hits))
        else:
            kegg_ids_list.append(None)
            annotations_list.append(None)
            pathways_list.append(None)
    
    result_df['kegg_ids'] = kegg_ids_list
    result_df['annotation'] = annotations_list
    result_df['pathway'] = pathways_list
    
    logger.info(f"Annotated {len(result_df[result_df['kegg_ids'].notna()])} genes with KEGG")
    
    return result_df

