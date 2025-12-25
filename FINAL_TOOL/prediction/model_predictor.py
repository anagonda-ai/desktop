"""
Model Predictor
Loads the combined trained XGBoost model and makes predictions for MGC candidates.
"""

import pandas as pd
import numpy as np
from typing import Dict, Optional
import os
import logging

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from FINAL_TOOL import config

try:
    import joblib
    JOBLIB_AVAILABLE = True
except ImportError:
    JOBLIB_AVAILABLE = False
    logger = logging.getLogger(__name__)
    logger.warning("joblib not available. Cannot load saved model. Install with: pip install joblib")

try:
    from xgboost import XGBClassifier
    XGBOOST_AVAILABLE = True
except ImportError:
    XGBOOST_AVAILABLE = False
    logger = logging.getLogger(__name__)
    logger.warning("xgboost not available. Cannot use XGBoost model. Install with: pip install xgboost")

logger = logging.getLogger(__name__)


class ModelPredictor:
    """
    Predictor for MGC candidates using the combined trained XGBoost model.
    """
    
    def __init__(self, weights_dir: str = None):
        """
        Initialize predictor with model directory.
        
        Args:
            weights_dir: Directory containing model files (default from config)
        """
        if weights_dir is None:
            weights_dir = config.MODEL_WEIGHTS_DIR
        
        self.weights_dir = weights_dir
        self.model = None
        self.feature_order = None
        
    def load_model(self) -> bool:
        """
        Load the saved XGBoost model and feature order.
        
        Returns:
            True if loaded successfully, False otherwise
        """
        if not JOBLIB_AVAILABLE:
            logger.error("joblib is required to load the model. Install with: pip install joblib")
            return False
        
        if not XGBOOST_AVAILABLE:
            logger.error("xgboost is required to use the model. Install with: pip install xgboost")
            return False
        
        model_file = os.path.join(self.weights_dir, config.COMBINED_MODEL)
        feature_order_file = os.path.join(self.weights_dir, config.COMBINED_FEATURE_ORDER)
        
        if not os.path.exists(model_file):
            logger.error(f"Model file not found: {model_file}")
            return False
        
        try:
            # Load the trained XGBoost model
            self.model = joblib.load(model_file)
            logger.info(f"Loaded XGBoost model from {model_file}")
            
            # Load feature order
            if os.path.exists(feature_order_file):
                with open(feature_order_file, 'r') as f:
                    self.feature_order = [line.strip() for line in f if line.strip()]
                logger.info(f"Loaded feature order ({len(self.feature_order)} features)")
            else:
                logger.warning(f"Feature order file not found: {feature_order_file}")
                # Try to get feature names from model
                if hasattr(self.model, 'feature_names_in_'):
                    self.feature_order = list(self.model.feature_names_in_)
                    logger.info(f"Using feature names from model ({len(self.feature_order)} features)")
                else:
                    logger.error("Cannot determine feature order")
                    return False
            
            return True
            
        except Exception as e:
            logger.error(f"Error loading model from {model_file}: {e}")
            return False
    
    def predict_mgc(self, features: Dict[str, float],
                   threshold: float = None) -> Dict[str, any]:
        """
        Predict MGC probability and classification using the trained XGBoost model.
        
        Args:
            features: Dictionary of feature values (with prefixes: promoter_, docking_, etc.)
            threshold: Classification threshold (default from config)
            
        Returns:
            Dictionary with combined probability and MGC prediction
        """
        if threshold is None:
            threshold = config.PREDICTION_THRESHOLD
        
        # Load model if not already loaded
        if self.model is None:
            if not self.load_model():
                return {
                    'combined_probability': np.nan,
                    'is_mgc': False,
                    'threshold': threshold,
                    'error': 'Failed to load model'
                }
        
        # Extract feature values in the order expected by the model
        feature_values = []
        missing_features = []
        
        for feature_name in self.feature_order:
            value = features.get(feature_name)
            
            if value is None or np.isnan(value):
                missing_features.append(feature_name)
                # XGBoost handles NaN/missing values natively
                value = np.nan
            
            feature_values.append(value)
        
        if missing_features:
            logger.warning(f"Missing {len(missing_features)} features: {missing_features[:5]}...")
        
        # Convert to numpy array (XGBoost handles NaN values)
        feature_array = np.array(feature_values).reshape(1, -1)
        
        # Make prediction using XGBoost model
        try:
            # XGBoost predict_proba returns probabilities for both classes [prob_class_0, prob_class_1]
            probabilities = self.model.predict_proba(feature_array)
            probability = float(probabilities[0][1])  # Probability of class 1 (MGC)
        except Exception as e:
            logger.error(f"Error making prediction: {e}")
            return {
                'combined_probability': np.nan,
                'is_mgc': False,
                'threshold': threshold,
                'error': str(e)
            }
        
        # Make classification
        is_mgc = probability >= threshold
        
        return {
            'combined_probability': probability,
            'is_mgc': bool(is_mgc),
            'threshold': threshold,
            'num_features_used': len(self.feature_order),
            'num_features_missing': len(missing_features)
        }
