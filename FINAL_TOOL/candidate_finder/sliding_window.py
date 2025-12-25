"""
Sliding Window Candidate Finder
Identifies potential MGC candidates using sliding windows over annotated genes.
"""

import pandas as pd
from typing import List, Dict, Tuple, Optional
from itertools import combinations, product
import logging

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from FINAL_TOOL import config

logger = logging.getLogger(__name__)


def clean_and_process_genes_and_pathways(genes_and_pathways: Dict[str, List[str]]) -> Dict[str, List[List[str]]]:
    """
    Clean and process pathway assignments for genes.
    Handles comma-separated pathways in pathway strings.
    
    Args:
        genes_and_pathways: Dictionary mapping gene_id to list of pathway strings
        
    Returns:
        Dictionary mapping gene_id to list of pathway lists
    """
    processed = {}
    for gene, pathways in genes_and_pathways.items():
        processed_paths = []
        for pathway_str in pathways:
            if pd.notna(pathway_str) and pathway_str:
                # Split by comma if multiple pathways
                pathway_list = [p.strip() for p in str(pathway_str).split(',') if p.strip()]
                if pathway_list:
                    processed_paths.append(pathway_list)
        if processed_paths:
            processed[gene] = processed_paths
    return processed


def find_first_common_element(genes_and_pathways: Dict[str, List[str]], min_genes: int) -> Tuple[Optional[str], List[str]]:
    """
    Find the first common pathway among genes that satisfies the minimum gene requirement.
    
    Args:
        genes_and_pathways: Dictionary mapping gene_id to list of pathway strings
        min_genes: Minimum number of genes required
        
    Returns:
        Tuple of (pathway, list of gene_ids) or (None, []) if no common pathway found
    """
    processed = clean_and_process_genes_and_pathways(genes_and_pathways)
    
    if len(processed) < min_genes:
        return None, []
    
    # Create sets of pathways for each gene
    gene_sets = {}
    for gene, paths in processed.items():
        # Flatten nested lists and create set
        flattened = set(item.strip() for sublist in paths for item in sublist)
        gene_sets[gene] = flattened
    
    # Iterate from largest possible group down to min_genes
    for group_size in range(len(gene_sets.keys()), min_genes - 1, -1):
        gene_combinations = combinations(gene_sets.keys(), group_size)
        for selected_genes in gene_combinations:
            selected_sets = [gene_sets[gene] for gene in selected_genes]
            
            # Find common pathway
            common_pathways = set.intersection(*selected_sets) if selected_sets else set()
            
            if common_pathways:
                # Return first common pathway and the gene list
                return list(common_pathways)[0], list(selected_genes)
    
    return None, []


def find_mgc_candidates(annotated_df: pd.DataFrame,
                        window_size: int = None,
                        min_metabolic_genes: int = None) -> pd.DataFrame:
    """
    Find MGC candidates using sliding window approach.
    
    Args:
        annotated_df: DataFrame with columns: id, chromosome, start, end, pathway, annotation, etc.
        window_size: Size of sliding window in genes (default from config)
        min_metabolic_genes: Minimum number of metabolic genes per window (default from config)
        
    Returns:
        DataFrame with candidate information
    """
    if window_size is None:
        window_size = config.WINDOW_SIZE
    if min_metabolic_genes is None:
        min_metabolic_genes = config.MIN_METABOLIC_GENES
    
    logger.info(f"Finding MGC candidates: window_size={window_size}, min_metabolic_genes={min_metabolic_genes}")
    
    # Filter genes with pathway annotations
    filtered_df = annotated_df[annotated_df['pathway'].notna()].copy()
    
    if filtered_df.empty:
        logger.warning("No genes with pathway annotations found")
        return pd.DataFrame()
    
    # Add index for window calculations
    filtered_df['index'] = filtered_df.index
    
    candidates = []
    chromosomes = filtered_df['chromosome'].unique()
    prev_matches = set()
    
    for chromosome in chromosomes:
        chromosome_data = filtered_df[filtered_df['chromosome'] == chromosome].copy()
        chromosome_data = chromosome_data.sort_values('start')  # Ensure sorted by position
        chromosome_data.reset_index(drop=True, inplace=True)
        chromosome_data['local_index'] = chromosome_data.index
        
        num_genes = len(chromosome_data)
        
        # Sliding window
        i = 0
        while i < num_genes:
            # Start window at position i
            window_start_idx = chromosome_data.iloc[i]['index']
            window = [chromosome_data.iloc[i]]
            
            # Extend window until we exceed window_size genes
            for j in range(i + 1, num_genes):
                gene_idx = chromosome_data.iloc[j]['index']
                # Check if this gene is within window_size distance (by index)
                if (gene_idx - window_start_idx) < window_size:
                    window.append(chromosome_data.iloc[j])
                else:
                    break
            
            # Check if window has enough genes
            if len(window) >= min_metabolic_genes:
                window_df = pd.DataFrame(window)
                
                # Create genes_and_pathways dictionary
                genes_and_pathways = {}
                genes_and_annotations = {}
                for _, row in window_df.iterrows():
                    genes_and_pathways[row['id']] = [row['pathway']] if pd.notna(row['pathway']) else []
                    genes_and_annotations[row['id']] = [row['annotation']] if pd.notna(row.get('annotation')) else []
                
                # Find common pathway
                pathway, metabolic_genes = find_first_common_element(genes_and_pathways, min_metabolic_genes)
                
                if pathway:
                    # Get annotations for metabolic genes
                    metabolic_annotations = []
                    for gene in metabolic_genes:
                        anns = genes_and_annotations.get(gene, [])
                        metabolic_annotations.append(anns[0] if anns else '')
                    
                    # Create candidate entry
                    candidate_key = tuple(sorted(metabolic_genes))
                    if candidate_key not in prev_matches:
                        prev_matches.add(candidate_key)
                        
                        candidate = {
                            'candidate_id': f"CANDIDATE_{len(candidates)+1:05d}",
                            'chromosome': chromosome,
                            'pathway': pathway,
                            'genes': ','.join(window_df['id'].tolist()),
                            'metabolic_genes': ','.join(metabolic_genes),
                            'metabolic_genes_annotations': ','.join(metabolic_annotations),
                            'start': window_df['start'].min(),
                            'end': window_df['end'].max(),
                            'num_genes': len(window_df),
                            'num_metabolic_genes': len(metabolic_genes)
                        }
                        candidates.append(candidate)
            
            # Move to next position
            i += 1
    
    if candidates:
        candidates_df = pd.DataFrame(candidates)
        logger.info(f"Found {len(candidates_df)} MGC candidates")
        return candidates_df
    else:
        logger.info("No MGC candidates found")
        return pd.DataFrame(columns=['candidate_id', 'chromosome', 'pathway', 'genes', 'metabolic_genes',
                                    'metabolic_genes_annotations', 'start', 'end', 'num_genes', 'num_metabolic_genes'])

