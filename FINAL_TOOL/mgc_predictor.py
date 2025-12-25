#!/usr/bin/env python3
"""
MGC Predictor - Main Pipeline
Complete tool for predicting Metabolic Gene Clusters from GFF3 and FASTA files.
"""

import argparse
import logging
import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from FINAL_TOOL.parsers.gff3_parser import extract_sequences_from_fasta
from FINAL_TOOL.annotation.kegg_annotator import annotate_genes_with_kegg
from FINAL_TOOL.candidate_finder.sliding_window import find_mgc_candidates
from FINAL_TOOL.feature_extraction import promoter_features
from FINAL_TOOL.feature_extraction import docking_features
from FINAL_TOOL.feature_extraction import cladepp_features
from FINAL_TOOL.feature_extraction import foldseek_features
from FINAL_TOOL.feature_extraction import e2p2_features
from FINAL_TOOL.prediction.model_predictor import ModelPredictor
from FINAL_TOOL.utils.sequence_utils import load_genome_fasta
from FINAL_TOOL import config

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def extract_all_features(candidate_genes: list,
                         genes_df: pd.DataFrame,
                         genome: dict,
                         work_dir: str = None) -> dict:
    """
    Extract all features for a candidate cluster.
    
    Args:
        candidate_genes: List of gene IDs in the candidate
        genes_df: DataFrame with gene information
        genome: Dictionary mapping chromosome to sequence
        work_dir: Working directory for temporary files
        
    Returns:
        Dictionary with all extracted features
    """
    all_features = {}
    
    # Extract promoter features
    logger.info("Extracting promoter features...")
    try:
        promoter_feats = promoter_features.extract_promoter_features(
            candidate_genes, genes_df, genome
        )
        # Add prefix
        for key, value in promoter_feats.items():
            all_features[f'promoter_{key}'] = value
    except Exception as e:
        logger.warning(f"Error extracting promoter features: {e}")
        for feat in ['mean_proximal_similarity', 'mean_distal_similarity', 
                    'mean_proximal_correlation', 'mean_distal_correlation', 'num_tfbs_types_found']:
            all_features[f'promoter_{feat}'] = np.nan
    
    # Extract docking features
    logger.info("Extracting docking features...")
    try:
        docking_feats = docking_features.extract_docking_features(
            candidate_genes, genes_df, work_dir
        )
        for key, value in docking_feats.items():
            all_features[f'docking_{key}'] = value
    except Exception as e:
        logger.warning(f"Error extracting docking features: {e}")
        for feat in ['mean_score_non_self', 'enrichment_score', 'z_score', 'effect_size',
                    'q75_score', 'fraction_strong_binders', 'fraction_weak_binders',
                    'max_score', 'median_score_non_self']:
            all_features[f'docking_{feat}'] = np.nan
    
    # Extract CladePP features
    logger.info("Extracting CladePP features...")
    try:
        cladepp_feats = cladepp_features.extract_cladepp_features(
            candidate_genes, genes_df, work_dir
        )
        for key, value in cladepp_feats.items():
            all_features[f'cladepp_{key}'] = value
    except Exception as e:
        logger.warning(f"Error extracting CladePP features: {e}")
        for feat in ['mean_cladepp_score', 'weighted_cladepp_score', 'positive_correlation_ratio',
                    'cladepp_multi_clade_high', 'cladepp_multi_clade_medium',
                    'cladepp_conservation_consistency', 'cladepp_max_pair_score']:
            all_features[f'cladepp_{feat}'] = np.nan
    
    # Extract Foldseek features
    logger.info("Extracting Foldseek features...")
    try:
        foldseek_feats = foldseek_features.extract_foldseek_features(
            candidate_genes, genes_df, work_dir
        )
        for key, value in foldseek_feats.items():
            all_features[f'foldseek_{key}'] = value
    except Exception as e:
        logger.warning(f"Error extracting Foldseek features: {e}")
        for feat in ['mean_score_non_self', 'enrichment_score', 'z_score', 'effect_size',
                    'foldseek_match_coverage']:
            all_features[f'foldseek_{feat}'] = np.nan
    
    # Extract E2P2 features
    logger.info("Extracting E2P2 features...")
    try:
        e2p2_feats = e2p2_features.extract_e2p2_features(
            candidate_genes, genes_df, work_dir
        )
        for key, value in e2p2_feats.items():
            all_features[f'e2p2_{key}'] = value
    except Exception as e:
        logger.warning(f"Error extracting E2P2 features: {e}")
        for feat in ['num_distinct_enzyme_classes', 'num_distinct_enzyme_subclasses',
                    'num_distinct_enzyme_families', 'num_distinct_enzyme_subfamilies',
                    'total_ec_numbers']:
            all_features[f'e2p2_{feat}'] = np.nan
    
    # Add cluster size
    all_features['cluster_size'] = len(candidate_genes)
    
    return all_features


def main():
    """Main pipeline function."""
    parser = argparse.ArgumentParser(
        description='Predict Metabolic Gene Clusters from GFF3 and FASTA files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python mgc_predictor.py genome.gff3 genome.fasta -o results.csv
  python mgc_predictor.py genome.gff3 genome.fasta -o results.csv --window-size 15 --min-genes 4
        """
    )
    
    parser.add_argument('gff3_file', type=str, help='Path to GFF3 annotation file')
    parser.add_argument('fasta_file', type=str, help='Path to genome FASTA file')
    parser.add_argument('-o', '--output', type=str, required=True,
                       help='Output CSV file path')
    parser.add_argument('--window-size', type=int, default=None,
                       help=f'Sliding window size (default: {config.WINDOW_SIZE})')
    parser.add_argument('--min-genes', type=int, default=None,
                       help=f'Minimum metabolic genes per window (default: {config.MIN_METABOLIC_GENES})')
    parser.add_argument('--temp-dir', type=str, default=None,
                       help=f'Temporary directory (default: {config.TEMP_DIR})')
    parser.add_argument('--kegg-db', type=str, default=None,
                       help=f'KEGG database path (default: {config.KEGG_DB})')
    parser.add_argument('--threshold', type=float, default=None,
                       help=f'Prediction threshold (default: {config.PREDICTION_THRESHOLD})')
    
    args = parser.parse_args()
    
    # Update config if provided
    if args.window_size:
        config.WINDOW_SIZE = args.window_size
    if args.min_genes:
        config.MIN_METABOLIC_GENES = args.min_genes
    if args.temp_dir:
        config.TEMP_DIR = args.temp_dir
    if args.kegg_db:
        config.KEGG_DB = args.kegg_db
    if args.threshold:
        config.PREDICTION_THRESHOLD = args.threshold
    
    logger.info("="*80)
    logger.info("MGC Predictor - Starting Pipeline")
    logger.info("="*80)
    logger.info(f"GFF3 file: {args.gff3_file}")
    logger.info(f"FASTA file: {args.fasta_file}")
    logger.info(f"Output file: {args.output}")
    logger.info(f"Window size: {config.WINDOW_SIZE}")
    logger.info(f"Min metabolic genes: {config.MIN_METABOLIC_GENES}")
    
    # Step 1: Parse GFF3 and extract sequences
    logger.info("\n" + "="*80)
    logger.info("Step 1: Parsing GFF3 and extracting sequences")
    logger.info("="*80)
    try:
        genes_df = extract_sequences_from_fasta(args.gff3_file, args.fasta_file)
        logger.info(f"✓ Parsed {len(genes_df)} genes")
    except Exception as e:
        logger.error(f"Error parsing GFF3: {e}")
        sys.exit(1)
    
    # Step 2: Annotate with KEGG
    logger.info("\n" + "="*80)
    logger.info("Step 2: Annotating genes with KEGG")
    logger.info("="*80)
    try:
        annotated_df = annotate_genes_with_kegg(genes_df)
        logger.info(f"✓ Annotated {len(annotated_df[annotated_df['kegg_ids'].notna()])} genes with KEGG")
    except Exception as e:
        logger.error(f"Error annotating with KEGG: {e}")
        sys.exit(1)
    
    # Step 3: Find MGC candidates
    logger.info("\n" + "="*80)
    logger.info("Step 3: Finding MGC candidates using sliding windows")
    logger.info("="*80)
    try:
        candidates_df = find_mgc_candidates(annotated_df)
        logger.info(f"✓ Found {len(candidates_df)} candidate clusters")
    except Exception as e:
        logger.error(f"Error finding candidates: {e}")
        sys.exit(1)
    
    if candidates_df.empty:
        logger.warning("No candidates found. Exiting.")
        candidates_df.to_csv(args.output, index=False)
        return
    
    # Step 4: Extract features and make predictions
    logger.info("\n" + "="*80)
    logger.info("Step 4: Extracting features and making predictions")
    logger.info("="*80)
    
    # Load genome for feature extraction
    genome = load_genome_fasta(args.fasta_file)
    
    # Initialize predictor
    predictor = ModelPredictor()
    
    # Process each candidate
    results = []
    for idx, candidate_row in candidates_df.iterrows():
        candidate_id = candidate_row['candidate_id']
        logger.info(f"\nProcessing candidate {idx+1}/{len(candidates_df)}: {candidate_id}")
        
        # Get gene list
        gene_list = candidate_row['metabolic_genes'].split(',')
        
        # Extract all features
        try:
            features = extract_all_features(gene_list, annotated_df, genome, config.TEMP_DIR)
        except Exception as e:
            logger.warning(f"Error extracting features for {candidate_id}: {e}")
            continue
        
        # Make prediction
        try:
            prediction = predictor.predict_mgc(features, threshold=config.PREDICTION_THRESHOLD)
        except Exception as e:
            logger.warning(f"Error predicting for {candidate_id}: {e}")
            continue
        
        # Create result row
        result_row = {
            'candidate_id': candidate_id,
            'chromosome': candidate_row['chromosome'],
            'pathway': candidate_row['pathway'],
            'genes': candidate_row['genes'],
            'metabolic_genes': candidate_row['metabolic_genes'],
            'start': candidate_row['start'],
            'end': candidate_row['end'],
            'num_genes': candidate_row['num_genes'],
            'num_metabolic_genes': candidate_row['num_metabolic_genes'],
            'is_mgc': prediction['is_mgc'],
            'combined_probability': prediction['combined_probability'],
            'threshold': prediction['threshold']
        }
        
        # Add prediction metadata
        if 'num_features_used' in prediction:
            result_row['num_features_used'] = prediction['num_features_used']
        if 'num_features_missing' in prediction:
            result_row['num_features_missing'] = prediction['num_features_missing']
        
        # Add all features
        result_row.update(features)
        
        results.append(result_row)
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # Step 5: Save results
    logger.info("\n" + "="*80)
    logger.info("Step 5: Saving results")
    logger.info("="*80)
    
    results_df.to_csv(args.output, index=False)
    logger.info(f"✓ Saved results to {args.output}")
    
    # Summary
    logger.info("\n" + "="*80)
    logger.info("Summary")
    logger.info("="*80)
    logger.info(f"Total candidates: {len(results_df)}")
    logger.info(f"Predicted MGCs: {len(results_df[results_df['is_mgc'] == True])}")
    logger.info(f"Predicted non-MGCs: {len(results_df[results_df['is_mgc'] == False])}")
    
    if len(results_df) > 0:
        logger.info(f"Mean probability: {results_df['combined_probability'].mean():.3f}")
        logger.info(f"Max probability: {results_df['combined_probability'].max():.3f}")
        logger.info(f"Min probability: {results_df['combined_probability'].min():.3f}")


if __name__ == "__main__":
    main()

