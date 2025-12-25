"""
Docking Feature Extraction
Extracts protein-protein docking features for MGC candidates from LightDock results.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Optional
import logging
import os
from pathlib import Path

logger = logging.getLogger(__name__)


def extract_docking_features_from_tsv(tsv_file: Path) -> Dict[str, float]:
    """
    Extract docking features from a LightDock all_vs_all.tsv results file.
    
    Args:
        tsv_file: Path to the TSV file containing docking results
        
    Returns:
        Dictionary with feature values
    """
    if not tsv_file.exists():
        logger.warning(f"Docking TSV file not found: {tsv_file}")
        return _get_default_docking_features()
    
    try:
        # Read the results file (tab-separated, no header)
        df = pd.read_csv(tsv_file, sep='\t', header=None)
        
        if len(df) == 0:
            logger.warning(f"Empty docking results file: {tsv_file}")
            return _get_default_docking_features()
        
        # Column 10 (index 10) is typically the score column in LightDock results
        # Adjust if needed based on actual format
        score_col = 10 if len(df.columns) > 10 else df.columns[-1]
        scores = df.iloc[:, score_col].astype(float)
        
        # Filter non-self comparisons (if first two columns are protein IDs)
        if len(df.columns) >= 2:
            non_self_mask = df.iloc[:, 0] != df.iloc[:, 1]
            scores_non_self = scores[non_self_mask] if non_self_mask.any() else scores
        else:
            scores_non_self = scores
        
        # Calculate basic statistics
        mean_score_non_self = float(scores_non_self.mean()) if len(scores_non_self) > 0 else 0.0
        median_score_non_self = float(scores_non_self.median()) if len(scores_non_self) > 0 else 0.0
        max_score = float(scores.max())
        q75_score = float(scores.quantile(0.75))
        
        # Calculate binding fractions based on score thresholds
        # LightDock scores are typically negative, where more negative = better binding
        # Thresholds: < -50 = strong, -50 to -25 = medium, >= -25 = weak
        strong_binders = (scores < -50).sum()
        medium_binders = ((scores >= -50) & (scores < -25)).sum()
        weak_binders = (scores >= -25).sum()
        total = len(scores)
        
        fraction_strong_binders = float(strong_binders / total) if total > 0 else 0.0
        fraction_weak_binders = float(weak_binders / total) if total > 0 else 0.0
        
        # Calculate enrichment and z-score (compared to random expectation)
        # Expected random score is typically around 0 or slightly negative
        expected_random = -10.0  # Conservative estimate for random interactions
        std_score = float(scores_non_self.std()) if len(scores_non_self) > 1 else 1.0
        
        enrichment_score = mean_score_non_self / expected_random if expected_random != 0 else 0.0
        z_score = (mean_score_non_self - expected_random) / std_score if std_score > 0 else 0.0
        effect_size = (mean_score_non_self - expected_random) / std_score if std_score > 0 else 0.0
        
        return {
            'mean_score_non_self': mean_score_non_self,
            'median_score_non_self': median_score_non_self,
            'max_score': max_score,
            'q75_score': q75_score,
            'fraction_strong_binders': fraction_strong_binders,
            'fraction_weak_binders': fraction_weak_binders,
            'enrichment_score': enrichment_score,
            'z_score': z_score,
            'effect_size': effect_size
        }
        
    except Exception as e:
        logger.error(f"Error parsing docking TSV file {tsv_file}: {e}")
        return _get_default_docking_features()


def extract_docking_features(candidate_genes: List[str],
                             genes_df: pd.DataFrame,
                             work_dir: Optional[str] = None,
                             docking_results_dir: Optional[str] = None) -> Dict[str, float]:
    """
    Extract docking features for a candidate cluster.
    
    This function looks for pre-computed LightDock results. If available,
    it parses them to compute features. Otherwise, returns NaN values.
    
    Args:
        candidate_genes: List of gene IDs in the candidate cluster
        genes_df: DataFrame with gene information
        work_dir: Working directory for temporary files
        docking_results_dir: Optional directory containing docking results
        
    Returns:
        Dictionary with feature values (NaN if results not available)
    """
    # If docking_results_dir is provided, look for results there
    if docking_results_dir:
        results_path = Path(docking_results_dir)
        
        # Try to find the all_vs_all.tsv file
        # Common patterns: results_tsv/{cluster_name}_all_vs_all.tsv
        # Or: {cluster_name}/results_tsv/{cluster_name}_all_vs_all.tsv
        
        # Try multiple possible locations
        cluster_name = "_".join(candidate_genes[:3]) if candidate_genes else "unknown"
        
        possible_paths = [
            results_path / "results_tsv" / f"{cluster_name}_all_vs_all.tsv",
            results_path / cluster_name / "results_tsv" / f"{cluster_name}_all_vs_all.tsv",
            results_path / f"{cluster_name}_all_vs_all.tsv",
        ]
        
        for tsv_path in possible_paths:
            if tsv_path.exists():
                logger.info(f"Found docking results at: {tsv_path}")
                return extract_docking_features_from_tsv(tsv_path)
    
    # No results found, return NaN values
    logger.warning("No docking results found. Docking features will be NaN.")
    return _get_default_docking_features()


def _get_default_docking_features() -> Dict[str, float]:
    """Return default docking features with NaN values"""
    return {
        'mean_score_non_self': np.nan,
        'median_score_non_self': np.nan,
        'max_score': np.nan,
        'q75_score': np.nan,
        'fraction_strong_binders': np.nan,
        'fraction_weak_binders': np.nan,
        'enrichment_score': np.nan,
        'z_score': np.nan,
        'effect_size': np.nan
    }
