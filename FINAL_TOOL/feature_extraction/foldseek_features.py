"""
Foldseek Feature Extraction
Extracts Foldseek structural matching features for MGC candidates from Foldseek results.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Optional
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def extract_foldseek_features_from_tsv(tsv_file: Path) -> Dict[str, float]:
    """
    Extract Foldseek features from a Foldseek all_vs_all.tsv results file.
    
    Args:
        tsv_file: Path to the TSV file containing Foldseek results
        
    Returns:
        Dictionary with feature values
    """
    if not tsv_file.exists():
        logger.warning(f"Foldseek TSV file not found: {tsv_file}")
        return _get_default_foldseek_features()
    
    try:
        # Read the results file (tab-separated, no header)
        df = pd.read_csv(tsv_file, sep='\t', header=None)
        
        if len(df) == 0:
            logger.warning(f"Empty Foldseek results file: {tsv_file}")
            return _get_default_foldseek_features()
        
        # Column 2 (index 2) is typically the TM-score in Foldseek results
        score_col = 2 if len(df.columns) > 2 else df.columns[-1]
        scores = df.iloc[:, score_col].astype(float)
        
        # Filter non-self comparisons
        if len(df.columns) >= 2:
            non_self_mask = df.iloc[:, 0] != df.iloc[:, 1]
            scores_non_self = scores[non_self_mask] if non_self_mask.any() else scores
        else:
            scores_non_self = scores
        
        # Calculate basic statistics
        mean_score_non_self = float(scores_non_self.mean()) if len(scores_non_self) > 0 else 0.0
        std_score_non_self = float(scores_non_self.std()) if len(scores_non_self) > 1 else 0.0
        
        # Calculate binding fractions based on TM-score thresholds
        # TM-score > 0.5 = strong structural match, > 0.3 = moderate, <= 0.3 = weak
        strong_binders = (scores > 0.5).sum()
        moderate_binders = ((scores > 0.3) & (scores <= 0.5)).sum()
        weak_binders = (scores <= 0.3).sum()
        total = len(scores)
        
        # Calculate match coverage: fraction of proteins with at least one strong match
        if len(df.columns) >= 2 and total > 0:
            strong_pairs = df[(df.iloc[:, 0] != df.iloc[:, 1]) & (df.iloc[:, score_col].astype(float) > 0.5)]
            if len(strong_pairs) > 0:
                proteins_with_strong_match = set(strong_pairs.iloc[:, 0]).union(set(strong_pairs.iloc[:, 1]))
                all_proteins = set(df.iloc[:, 0]).union(set(df.iloc[:, 1]))
                foldseek_match_coverage = len(proteins_with_strong_match) / len(all_proteins) if len(all_proteins) > 0 else 0.0
            else:
                foldseek_match_coverage = 0.0
        else:
            foldseek_match_coverage = 0.0
        
        # Calculate enrichment and statistical measures
        expected_random = 0.3  # Expected random TM-score
        enrichment_score = mean_score_non_self / expected_random if expected_random > 0 else 0.0
        z_score = (mean_score_non_self - expected_random) / std_score_non_self if std_score_non_self > 0 else 0.0
        effect_size = (mean_score_non_self - expected_random) / std_score_non_self if std_score_non_self > 0 else 0.0
        
        return {
            'mean_score_non_self': mean_score_non_self,
            'enrichment_score': enrichment_score,
            'z_score': z_score,
            'effect_size': effect_size,
            'foldseek_match_coverage': foldseek_match_coverage
        }
        
    except Exception as e:
        logger.error(f"Error parsing Foldseek TSV file {tsv_file}: {e}")
        return _get_default_foldseek_features()


def extract_foldseek_features(candidate_genes: List[str],
                              genes_df: pd.DataFrame,
                              work_dir: Optional[str] = None,
                              foldseek_results_dir: Optional[str] = None) -> Dict[str, float]:
    """
    Extract Foldseek structural matching features for a candidate cluster.
    
    This function looks for pre-computed Foldseek results. If available,
    it parses them to compute features. Otherwise, returns NaN values.
    
    Args:
        candidate_genes: List of gene IDs in the candidate cluster
        genes_df: DataFrame with gene information
        work_dir: Working directory for temporary files
        foldseek_results_dir: Optional directory containing Foldseek results
        
    Returns:
        Dictionary with feature values (NaN if results not available)
    """
    # If foldseek_results_dir is provided, look for results there
    if foldseek_results_dir:
        results_path = Path(foldseek_results_dir)
        
        # Try to find the all_vs_all.tsv file
        cluster_name = "_".join(candidate_genes[:3]) if candidate_genes else "unknown"
        
        possible_paths = [
            results_path / f"{cluster_name}_all_vs_all.tsv",
            results_path / cluster_name / f"{cluster_name}_all_vs_all.tsv",
        ]
        
        for tsv_path in possible_paths:
            if tsv_path.exists():
                logger.info(f"Found Foldseek results at: {tsv_path}")
                return extract_foldseek_features_from_tsv(tsv_path)
    
    # No results found, return NaN values
    logger.warning("No Foldseek results found. Foldseek features will be NaN.")
    return _get_default_foldseek_features()


def _get_default_foldseek_features() -> Dict[str, float]:
    """Return default Foldseek features with NaN values"""
    return {
        'mean_score_non_self': np.nan,
        'enrichment_score': np.nan,
        'z_score': np.nan,
        'effect_size': np.nan,
        'foldseek_match_coverage': np.nan
    }
