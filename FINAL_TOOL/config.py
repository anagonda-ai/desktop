"""
Configuration file for MGC Prediction Tool
"""

import os
from pathlib import Path

# Get FINAL_TOOL directory (parent of this config.py file)
FINAL_TOOL_DIR = Path(__file__).parent

# Paths
KEGG_DB = str(FINAL_TOOL_DIR / "data" / "merged_metabolic_pathways.fasta")
MODEL_WEIGHTS_DIR = str(FINAL_TOOL_DIR / "models")

# Combined model files (relative to MODEL_WEIGHTS_DIR)
# This is the final trained XGBoost model with all features combined
COMBINED_MODEL = "combined_model.pkl"  # Full XGBoost model
COMBINED_WEIGHTS = "combined_all_features_weights.csv"  # Feature importances (for reference)
COMBINED_FEATURE_ORDER = "combined_model_feature_order.txt"  # Feature order for prediction

# Sliding window parameters
WINDOW_SIZE = 10
MIN_METABOLIC_GENES = 3

# BLAST parameters
BLAST_EVALUE = 1e-3
BLAST_IDENTITY = 50.0
BLAST_COVERAGE = 70.0

# Promoter extraction
PROMOTER_UPSTREAM = 1000

# Prediction threshold (F1-optimal from training)
# These will be loaded from model summary files if available
PREDICTION_THRESHOLD = 0.5

# Temporary directory for intermediate files
TEMP_DIR = "/tmp/mgc_predictor"

