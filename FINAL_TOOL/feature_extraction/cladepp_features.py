"""
CladePP Feature Extraction
Extracts CladePP phylogenetic profiling features for MGC candidates from CladePP results.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Optional
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def extract_cladepp_features_from_summary(summary_file: Path) -> Dict[str, float]:
    """
    Extract CladePP features from a summary.csv file.
    
    Args:
        summary_file: Path to the CladePP summary.csv file
        
    Returns:
        Dictionary with feature values
    """
    if not summary_file.exists():
        logger.warning(f"CladePP summary file not found: {summary_file}")
        return _get_default_cladepp_features()
    
    try:
        # Read the summary CSV
        df = pd.read_csv(summary_file)
        
        if len(df) == 0:
            logger.warning(f"Empty CladePP summary file: {summary_file}")
            return _get_default_cladepp_features()
        
        # Filter successful analyses
        if 'status' in df.columns:
            df = df[df['status'] == 'success'].copy()
        
        if len(df) == 0:
            logger.warning(f"No successful CladePP analyses in: {summary_file}")
            return _get_default_cladepp_features()
        
        # Get cladepp scores
        if 'cladepp_score' not in df.columns:
            logger.warning(f"'cladepp_score' column not found in: {summary_file}")
            return _get_default_cladepp_features()
        
        scores = df['cladepp_score'].astype(float).dropna()
        
        if len(scores) == 0:
            logger.warning(f"No valid cladepp_score values in: {summary_file}")
            return _get_default_cladepp_features()
        
        # Calculate mean cladepp score
        mean_cladepp_score = float(scores.mean())
        
        # Calculate weighted cladepp score (weighted by clade size if available)
        if 'clade_size' in df.columns:
            sizes = df['clade_size'].astype(float)
            weights = sizes / sizes.sum() if sizes.sum() > 0 else np.ones(len(scores)) / len(scores)
            weighted_cladepp_score = float((scores * weights[:len(scores)]).sum())
        else:
            weighted_cladepp_score = mean_cladepp_score
        
        # Calculate positive correlation ratio
        if 'positive_corr_pairs' in df.columns and 'total_corr_pairs' in df.columns:
            positive_pairs = df['positive_corr_pairs'].fillna(0).sum()
            total_pairs = df['total_corr_pairs'].fillna(0).sum()
            positive_correlation_ratio = positive_pairs / total_pairs if total_pairs > 0 else 0.0
        elif 'mean_anchor_corr' in df.columns:
            # Use mean anchor correlation as proxy
            mean_corr = df['mean_anchor_corr'].astype(float).dropna()
            positive_correlation_ratio = float((mean_corr > 0).mean()) if len(mean_corr) > 0 else 0.0
        else:
            positive_correlation_ratio = 0.0
        
        # Calculate multi-clade conservation metrics
        # High threshold: score > 0.8, medium threshold: score > 0.6
        cladepp_multi_clade_high = float((scores > 0.8).mean())
        cladepp_multi_clade_medium = float((scores > 0.6).mean())
        
        # Calculate conservation consistency (inverse of coefficient of variation)
        mean_score = scores.mean()
        std_score = scores.std()
        if mean_score > 0:
            consistency = 1 - (std_score / mean_score)
            cladepp_conservation_consistency = max(0.0, min(1.0, float(consistency)))
        else:
            cladepp_conservation_consistency = 0.0
        
        # Maximum pair co-evolution score
        cladepp_max_pair_score = float(scores.max())
        
        return {
            'mean_cladepp_score': mean_cladepp_score,
            'weighted_cladepp_score': weighted_cladepp_score,
            'positive_correlation_ratio': positive_correlation_ratio,
            'cladepp_multi_clade_high': cladepp_multi_clade_high,
            'cladepp_multi_clade_medium': cladepp_multi_clade_medium,
            'cladepp_conservation_consistency': cladepp_conservation_consistency,
            'cladepp_max_pair_score': cladepp_max_pair_score
        }
        
    except Exception as e:
        logger.error(f"Error parsing CladePP summary file {summary_file}: {e}")
        return _get_default_cladepp_features()


def extract_cladepp_features(candidate_genes: List[str],
                             genes_df: pd.DataFrame,
                             work_dir: Optional[str] = None,
                             cladepp_results_dir: Optional[str] = None) -> Dict[str, float]:
    """
    Extract CladePP phylogenetic profiling features for a candidate cluster.
    
    This function looks for pre-computed CladePP results. If available,
    it parses them to compute features. Otherwise, returns NaN values.
    
    Args:
        candidate_genes: List of gene IDs in the candidate cluster
        genes_df: DataFrame with gene information
        work_dir: Working directory for temporary files
        cladepp_results_dir: Optional directory containing CladePP results
        
    Returns:
        Dictionary with feature values (NaN if results not available)
    """
    # If cladepp_results_dir is provided, look for results there
    if cladepp_results_dir:
        results_path = Path(cladepp_results_dir)
        
        # Try to find the summary.csv file
        cluster_name = "_".join(candidate_genes[:3]) if candidate_genes else "unknown"
        
        possible_paths = [
            results_path / cluster_name / "summary.csv",
            results_path / "summary.csv",
            results_path / f"{cluster_name}_summary.csv",
        ]
        
        for summary_path in possible_paths:
            if summary_path.exists():
                logger.info(f"Found CladePP results at: {summary_path}")
                return extract_cladepp_features_from_summary(summary_path)
    
    # No results found, return NaN values
    logger.warning("No CladePP results found. CladePP features will be NaN.")
    return _get_default_cladepp_features()


def _get_default_cladepp_features() -> Dict[str, float]:
    """Return default CladePP features with NaN values"""
    return {
        'mean_cladepp_score': np.nan,
        'weighted_cladepp_score': np.nan,
        'positive_correlation_ratio': np.nan,
        'cladepp_multi_clade_high': np.nan,
        'cladepp_multi_clade_medium': np.nan,
        'cladepp_conservation_consistency': np.nan,
        'cladepp_max_pair_score': np.nan
    }
