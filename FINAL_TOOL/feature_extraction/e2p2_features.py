"""
E2P2 Feature Extraction
Extracts E2P2 enzyme classification features for MGC candidates from E2P2 results.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Optional
import logging
import re
from pathlib import Path

logger = logging.getLogger(__name__)


def parse_e2p2_default_pf(e2p2_file: Path) -> List[Dict]:
    """
    Parse E2P2 results from .default.pf file and extract EC numbers.
    
    Args:
        e2p2_file: Path to the .default.pf file
        
    Returns:
        List of dictionaries containing protein EC number information
    """
    results = []
    current_entry = {}
    
    if not e2p2_file.exists():
        logger.warning(f"E2P2 file not found: {e2p2_file}")
        return results
    
    try:
        with open(e2p2_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('ID\t'):
                    if current_entry and 'protein_id' in current_entry:
                        results.append(current_entry)
                    current_entry = {'protein_id': line[3:].strip()}
                elif line.startswith('EC\t'):
                    # Handle multiple EC numbers per protein
                    if 'ec_numbers' not in current_entry:
                        current_entry['ec_numbers'] = []
                    current_entry['ec_numbers'].append(line[3:].strip())
                elif line.startswith('NAME\t'):
                    current_entry['protein_name'] = line[5:].strip()
                elif line.startswith('PRODUCT-TYPE\t'):
                    current_entry['product_type'] = line[13:].strip()
                elif line == '//':
                    if current_entry and 'protein_id' in current_entry:
                        results.append(current_entry)
                    current_entry = {}
        
        # Add last entry if exists
        if current_entry and 'protein_id' in current_entry:
            results.append(current_entry)
            
    except Exception as e:
        logger.error(f"Error parsing E2P2 file {e2p2_file}: {e}")
    
    return results


def extract_ec_features_from_results(e2p2_results: List[Dict]) -> Dict[str, float]:
    """
    Extract EC number features from parsed E2P2 results.
    
    Args:
        e2p2_results: List of dictionaries from parse_e2p2_default_pf
        
    Returns:
        Dictionary with EC feature values
    """
    all_ec_numbers = []
    
    for result in e2p2_results:
        ec_numbers = result.get('ec_numbers', [])
        for ec_number in ec_numbers:
            if ec_number and ec_number != 'UNKNOWN':
                # Clean EC number (remove any extra text)
                ec_clean = re.sub(r'[^\d\.]', '', ec_number)
                if re.match(r'^\d+\.\d+\.\d+\.\d+$', ec_clean):
                    all_ec_numbers.append(ec_clean)
    
    # Extract enzyme classes at different levels
    enzyme_classes = set()          # First level: X
    enzyme_subclasses = set()       # First two levels: X.Y
    enzyme_families = set()         # First three levels: X.Y.Z
    enzyme_subfamilies = set()      # All levels: X.Y.Z.W
    
    for ec in all_ec_numbers:
        parts = ec.split('.')
        if len(parts) >= 1:
            enzyme_classes.add(parts[0])
        if len(parts) >= 2:
            enzyme_subclasses.add(f"{parts[0]}.{parts[1]}")
        if len(parts) >= 3:
            enzyme_families.add(f"{parts[0]}.{parts[1]}.{parts[2]}")
        if len(parts) >= 4:
            enzyme_subfamilies.add(f"{parts[0]}.{parts[1]}.{parts[2]}.{parts[3]}")
    
    return {
        'num_distinct_enzyme_classes': len(enzyme_classes),
        'num_distinct_enzyme_subclasses': len(enzyme_subclasses),
        'num_distinct_enzyme_families': len(enzyme_families),
        'num_distinct_enzyme_subfamilies': len(enzyme_subfamilies),
        'total_ec_numbers': len(all_ec_numbers)
    }


def extract_e2p2_features_from_file(e2p2_file: Path) -> Dict[str, float]:
    """
    Extract E2P2 features from a .default.pf file.
    
    Args:
        e2p2_file: Path to the E2P2 .default.pf file
        
    Returns:
        Dictionary with feature values
    """
    if not e2p2_file.exists():
        logger.warning(f"E2P2 file not found: {e2p2_file}")
        return _get_default_e2p2_features()
    
    try:
        # Parse E2P2 results
        e2p2_results = parse_e2p2_default_pf(e2p2_file)
        
        if len(e2p2_results) == 0:
            logger.warning(f"No E2P2 results in: {e2p2_file}")
            return _get_default_e2p2_features()
        
        # Extract EC features
        ec_features = extract_ec_features_from_results(e2p2_results)
        
        return ec_features
        
    except Exception as e:
        logger.error(f"Error extracting E2P2 features from {e2p2_file}: {e}")
        return _get_default_e2p2_features()


def extract_e2p2_features(candidate_genes: List[str],
                          genes_df: pd.DataFrame,
                          work_dir: Optional[str] = None,
                          e2p2_results_dir: Optional[str] = None) -> Dict[str, float]:
    """
    Extract E2P2 enzyme classification features for a candidate cluster.
    
    This function looks for pre-computed E2P2 results. If available,
    it parses them to compute features. Otherwise, returns NaN values.
    
    Args:
        candidate_genes: List of gene IDs in the candidate cluster
        genes_df: DataFrame with gene information
        work_dir: Working directory for temporary files
        e2p2_results_dir: Optional directory containing E2P2 results
        
    Returns:
        Dictionary with feature values (NaN if results not available)
    """
    # If e2p2_results_dir is provided, look for results there
    if e2p2_results_dir:
        results_path = Path(e2p2_results_dir)
        
        # Try to find the .default.pf file
        cluster_name = "_".join(candidate_genes[:3]) if candidate_genes else "unknown"
        
        possible_paths = [
            results_path / cluster_name / f"{cluster_name}.MaxWeightAbsoluteThreshold.default.pf",
            results_path / f"{cluster_name}.MaxWeightAbsoluteThreshold.default.pf",
            results_path / cluster_name / f"{cluster_name}.default.pf",
            results_path / f"{cluster_name}.default.pf",
        ]
        
        for e2p2_path in possible_paths:
            if e2p2_path.exists():
                logger.info(f"Found E2P2 results at: {e2p2_path}")
                return extract_e2p2_features_from_file(e2p2_path)
    
    # No results found, return NaN values
    logger.warning("No E2P2 results found. E2P2 features will be NaN.")
    return _get_default_e2p2_features()


def _get_default_e2p2_features() -> Dict[str, float]:
    """Return default E2P2 features with NaN values"""
    return {
        'num_distinct_enzyme_classes': np.nan,
        'num_distinct_enzyme_subclasses': np.nan,
        'num_distinct_enzyme_families': np.nan,
        'num_distinct_enzyme_subfamilies': np.nan,
        'total_ec_numbers': np.nan
    }
