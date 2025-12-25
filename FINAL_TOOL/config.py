"""
Configuration file for MGC Prediction Tool
"""

# Paths
KEGG_DB = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/fasta/merged_metabolic_pathways.fasta"
MODEL_WEIGHTS_DIR = "/groups/itay_mayrose/alongonda/desktop/mibig_validate"

# Combined model files (relative to MODEL_WEIGHTS_DIR)
# This is the final trained XGBoost model with all features combined
COMBINED_MODEL = "combined/combined_model.pkl"  # Full XGBoost model
COMBINED_WEIGHTS = "combined/combined_all_features_weights.csv"  # Feature importances (for reference)
COMBINED_FEATURE_ORDER = "combined/combined_model_feature_order.txt"  # Feature order for prediction

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

