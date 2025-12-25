"""
Promoter Feature Extraction
Extracts promoter similarity features for MGC candidates based on TFBS pattern analysis.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Optional
import logging
import re
from collections import defaultdict
from itertools import combinations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from FINAL_TOOL.utils.sequence_utils import load_fasta_genome, extract_promoter_sequence
from FINAL_TOOL import config

logger = logging.getLogger(__name__)


class PromoterFeatureExtractor:
    """
    Extracts promoter features using TFBS pattern analysis.
    Based on the implementation in python_scripts/features/extract_promotor/promoter_analysis.py
    """
    
    def __init__(self, upstream_length: int = config.PROMOTER_UPSTREAM):
        self.upstream_length = upstream_length
        
        # Plant-specific TFBS patterns from PlantCARE, PLACE, JASPAR databases
        # Same patterns as used in the original implementation
        self.plant_tfbs_patterns = {
            # Core promoter elements
            'TATA_box': [r'TATAAA', r'TATAWAW', r'TATAWAR'],
            'CAAT_box': [r'CAAT', r'CCAAT', r'CAAAT'],
            'Initiator': [r'[CT][CT]A[ATGC][AT][CT][CT]'],
            
            # Major plant TF families
            'MYB_binding': [r'[CT]AACNG', r'YAACKG', r'CAACAG'],
            'bZIP_binding': [r'[TG]GACGT[CA]', r'ACGT', r'CACGTG'],
            'WRKY_binding': [r'[CT]TGAC[CT]', r'TTGACY'],
            'AP2_ERF': [r'GCCGCC', r'AGCCGCC'],
            'bHLH_binding': [r'CANNTG', r'CACGTG', r'CATGTG'],
            'NAC_binding': [r'CACG', r'CATGTG'],
            'GATA_binding': [r'[AT]GATA[AG]', r'HGATAR'],
            'DOF_binding': [r'AAAGG', r'AAAG'],
            'TCP_binding': [r'GGNCCCAC'],
            'HSF_binding': [r'GAANNTTC', r'GAANNTTCNNGAANNTTC'],
            
            # Hormone response elements
            'ABA_response': [r'ACGTGG[CT]', r'ACGTGTC'],
            'Auxin_response': [r'TGTCTC', r'GAGACA'],
            'GA_response': [r'TAACAAA', r'TATCCAC'],
            'Cytokinin_response': [r'TATTAG', r'AGATCCT'],
            'Ethylene_response': [r'ATTTCAAA', r'TAAGAGCCGCC'],
            
            # Stress response elements
            'Drought_response': [r'TACCGACAT', r'ACGTGG[CT]'],
            'Cold_response': [r'TGGCCGAC', r'CCGAC'],
            'Heat_response': [r'GAANNTTC'],
            'Light_response': [r'ACGT', r'CAATCA', r'ATCTTA'],
            
            # Tissue-specific elements
            'Root_specific': [r'ATATT', r'TTATT'],
            'Leaf_specific': [r'CAATCA', r'ATCTTA'],
            'Seed_specific': [r'ACAAAA', r'CATGCA'],
            
            # Other regulatory elements
            'CpG_island': [r'CG'],
            'GC_box': [r'GGGCGG', r'CCGCCC'],
            'E_box': [r'CANNTG'],
            'CCAAT_box': [r'CCAAT'],
            'Silencer': [r'AATTT'],
            'Enhancer': [r'TGACGT[CA]']
        }
    
    def scan_tfbs_patterns(self, sequence: str) -> Dict[str, List[Dict]]:
        """
        Scan sequence for known TFBS patterns.
        
        Args:
            sequence: DNA sequence to scan
            
        Returns:
            Dictionary mapping TFBS type to list of match positions
        """
        sequence = str(sequence).upper()
        tfbs_matches = defaultdict(list)
        
        for tfbs_name, patterns in self.plant_tfbs_patterns.items():
            for pattern in patterns:
                # Convert IUPAC nucleotide codes to regex
                regex_pattern = pattern
                iupac_codes = {
                    'W': '[AT]', 'S': '[GC]', 'M': '[AC]', 'K': '[GT]',
                    'R': '[AG]', 'Y': '[CT]', 'B': '[CGT]', 'D': '[AGT]',
                    'H': '[ACT]', 'V': '[ACG]', 'N': '[ATGC]'
                }
                
                for code, replacement in iupac_codes.items():
                    regex_pattern = regex_pattern.replace(code, replacement)
                
                # Find all matches
                try:
                    matches = list(re.finditer(regex_pattern, sequence))
                    for match in matches:
                        tfbs_matches[tfbs_name].append({
                            'start': match.start(),
                            'end': match.end(),
                            'sequence': match.group()
                        })
                except re.error:
                    logger.warning(f"Invalid regex pattern for {tfbs_name}: {pattern}")
                    continue
        
        return dict(tfbs_matches)
    
    def calculate_tfbs_similarity(self, tfbs1: Dict[str, List], tfbs2: Dict[str, List]) -> tuple:
        """
        Calculate similarity between two TFBS profiles.
        
        Returns:
            Tuple of (jaccard_similarity, correlation)
        """
        # Get all TFBS types present in either promoter
        all_types = set(tfbs1.keys()) | set(tfbs2.keys())
        if not all_types:
            return 0.0, 0.0
        
        # Jaccard similarity: shared TFBS types / total TFBS types
        shared_types = set(tfbs1.keys()) & set(tfbs2.keys())
        jaccard_similarity = len(shared_types) / len(all_types) if all_types else 0.0
        
        # Correlation of TFBS densities
        counts1 = [len(tfbs1.get(tfbs_type, [])) for tfbs_type in all_types]
        counts2 = [len(tfbs2.get(tfbs_type, [])) for tfbs_type in all_types]
        
        if sum(counts1) > 0 and sum(counts2) > 0:
            # Normalize counts to get densities
            total1, total2 = sum(counts1), sum(counts2)
            norm_counts1 = np.array([c/total1 for c in counts1])
            norm_counts2 = np.array([c/total2 for c in counts2])
            
            # Calculate Pearson correlation
            if np.std(norm_counts1) > 0 and np.std(norm_counts2) > 0:
                correlation = np.corrcoef(norm_counts1, norm_counts2)[0, 1]
                correlation = max(0.0, float(correlation))  # Only positive correlations are meaningful
            else:
                correlation = 0.0
        else:
            correlation = 0.0
        
        return jaccard_similarity, correlation
    
    def analyze_promoter_regions(self, promoter_sequences: Dict[str, str]) -> Dict[str, Dict]:
        """
        Analyze proximal and distal promoter regions separately.
        
        Args:
            promoter_sequences: Dictionary mapping gene_id to promoter sequence
            
        Returns:
            Dictionary mapping gene_id to regional analysis data
        """
        regional_analysis = {}
        
        for gene_id, seq_str in promoter_sequences.items():
            seq_str = str(seq_str).upper()
            length = len(seq_str)
            
            # Split promoter into regions (TSS is at the end of the sequence for + strand)
            # Standard split: distal (-800 to -200), proximal (-200 to TSS)
            if length >= 800:
                # Standard split: distal (-800 to -200), proximal (-200 to TSS)
                distal_region = seq_str[:600]    # First 600bp = distal (-800 to -200)
                proximal_region = seq_str[600:]  # Last 200bp = proximal (-200 to TSS)
            elif length >= 400:
                # Shorter promoters: split in half
                split_point = length // 2
                distal_region = seq_str[:split_point]
                proximal_region = seq_str[split_point:]
            else:
                # Very short: use entire sequence for both
                distal_region = seq_str
                proximal_region = seq_str
            
            # Scan for TFBS in each region
            distal_tfbs = self.scan_tfbs_patterns(distal_region)
            proximal_tfbs = self.scan_tfbs_patterns(proximal_region)
            
            regional_analysis[gene_id] = {
                'total_length': length,
                'distal_region': {
                    'length': len(distal_region),
                    'tfbs': distal_tfbs,
                    'tfbs_count': sum(len(matches) for matches in distal_tfbs.values())
                },
                'proximal_region': {
                    'length': len(proximal_region),
                    'tfbs': proximal_tfbs,
                    'tfbs_count': sum(len(matches) for matches in proximal_tfbs.values())
                }
            }
        
        return regional_analysis
    
    def _get_default_features(self) -> Dict[str, float]:
        """Return default promoter features with NaN values"""
        return {
            'promoter_mean_proximal_similarity': np.nan,
            'promoter_mean_distal_similarity': np.nan,
            'promoter_mean_proximal_correlation': np.nan,
            'promoter_mean_distal_correlation': np.nan,
            'promoter_num_tfbs_types_found': 0
        }


# Main function for compatibility with existing code
def extract_promoter_features(candidate_genes: List[str],
                              genes_df: pd.DataFrame,
                              genome: Dict = None,
                              genome_fasta_path: str = None,
                              upstream: int = None) -> Dict[str, float]:
    """
    Extract promoter similarity features for a candidate cluster.
    
    Args:
        candidate_genes: List of gene IDs in the candidate cluster
        genes_df: DataFrame with gene information (must include gene_id, chromosome, start, end, strand)
        genome: Dictionary mapping chromosome to sequence (if already loaded)
        genome_fasta_path: Path to genome FASTA file (if genome not provided)
        upstream: Number of bp upstream to extract (default from config)
        
    Returns:
        Dictionary with feature values (without 'promoter_' prefix, as it will be added by caller)
    """
    if upstream is None:
        upstream = config.PROMOTER_UPSTREAM
    
    # Filter genes_df to only candidate genes
    candidate_genes_df = genes_df[genes_df['gene_id'].isin(candidate_genes)].copy()
    
    if candidate_genes_df.empty:
        logger.warning(f"No genes found for candidate cluster")
        extractor = PromoterFeatureExtractor(upstream_length=upstream)
        default_features = extractor._get_default_features()
        # Remove 'promoter_' prefix for return
        return {k.replace('promoter_', '') if k.startswith('promoter_') else k: v 
                for k, v in default_features.items()}
    
    # Determine genome source
    if genome is None and genome_fasta_path is None:
        logger.error("Either genome dict or genome_fasta_path must be provided")
        extractor = PromoterFeatureExtractor(upstream_length=upstream)
        default_features = extractor._get_default_features()
        return {k.replace('promoter_', '') if k.startswith('promoter_') else k: v 
                for k, v in default_features.items()}
    
    # Create extractor
    extractor = PromoterFeatureExtractor(upstream_length=upstream)
    
    # If genome dict provided, use it directly; otherwise load from file
    if genome is None:
        genome = load_fasta_genome(genome_fasta_path)
        if not genome:
            logger.error("Failed to load genome")
            default_features = extractor._get_default_features()
            return {k.replace('promoter_', '') if k.startswith('promoter_') else k: v 
                    for k, v in default_features.items()}
    
    # Extract promoter sequences
    promoter_sequences = {}
    for _, row in candidate_genes_df.iterrows():
        gene_id = row['gene_id']
        chromosome = row['chromosome']
        start = row['start']
        end = row['end']
        strand = row['strand']
        
        # Extract upstream region using the utility function
        promoter_seq_str = extract_promoter_sequence(genome, chromosome, start, end, strand, upstream=upstream)
        
        if promoter_seq_str:
            promoter_sequences[gene_id] = promoter_seq_str
    
    if len(promoter_sequences) < 2:
        logger.warning("Less than 2 promoter sequences extracted, cannot compute pairwise similarities.")
        default_features = extractor._get_default_features()
        return {k.replace('promoter_', '') if k.startswith('promoter_') else k: v 
                for k, v in default_features.items()}
    
    # Analyze promoter regions
    regional_analysis = extractor.analyze_promoter_regions(promoter_sequences)
    
    # Calculate pairwise similarities
    gene_ids = list(regional_analysis.keys())
    proximal_similarities = []
    distal_similarities = []
    proximal_correlations = []
    distal_correlations = []
    
    for i, j in combinations(range(len(gene_ids)), 2):
        gene1_id = gene_ids[i]
        gene2_id = gene_ids[j]
        
        res1 = regional_analysis[gene1_id]
        res2 = regional_analysis[gene2_id]
        
        # Proximal region
        jaccard_sim_prox, pearson_corr_prox = extractor.calculate_tfbs_similarity(
            res1['proximal_region']['tfbs'],
            res2['proximal_region']['tfbs']
        )
        proximal_similarities.append(jaccard_sim_prox)
        proximal_correlations.append(pearson_corr_prox)
        
        # Distal region
        jaccard_sim_dist, pearson_corr_dist = extractor.calculate_tfbs_similarity(
            res1['distal_region']['tfbs'],
            res2['distal_region']['tfbs']
        )
        distal_similarities.append(jaccard_sim_dist)
        distal_correlations.append(pearson_corr_dist)
    
    # Collect all TFBS types found
    all_tfbs_types = set()
    for data in regional_analysis.values():
        all_tfbs_types.update(data['proximal_region']['tfbs'].keys())
        all_tfbs_types.update(data['distal_region']['tfbs'].keys())
    
    features = {
        'mean_proximal_similarity': float(np.mean(proximal_similarities)) if proximal_similarities else np.nan,
        'mean_distal_similarity': float(np.mean(distal_similarities)) if distal_similarities else np.nan,
        'mean_proximal_correlation': float(np.mean(proximal_correlations)) if proximal_correlations else np.nan,
        'mean_distal_correlation': float(np.mean(distal_correlations)) if distal_correlations else np.nan,
        'num_tfbs_types_found': len(all_tfbs_types) if all_tfbs_types else 0
    }
    
    return features
