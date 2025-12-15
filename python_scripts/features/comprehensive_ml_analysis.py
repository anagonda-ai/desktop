#!/usr/bin/env python3
"""
Comprehensive ML Analysis - All 5 Feature Categories
Runs all ML classification models and creates comprehensive analysis
"""

import subprocess
import sys
import os
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (roc_curve, auc, precision_recall_curve, confusion_matrix, 
                            matthews_corrcoef, f1_score, make_scorer)
from sklearn.model_selection import (cross_validate, 
                                    learning_curve, validation_curve)
import warnings
import time
from scipy.stats import loguniform

# Ray framework for distributed ML - Required, no fallbacks
import ray
from ray import tune, train
from ray.tune.search.optuna import OptunaSearch
from ray.tune.schedulers import ASHAScheduler
from ray.data import from_pandas

warnings.filterwarnings('ignore')

# Data file paths for each category
DATA_FILES = {
    'CladePP': {
        'path': '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/cladepp_cluster_metrics_enhanced.csv',
        'cluster_col': 'name',
        'label_col': 'name'  # Will create label from this
    },
    'Docking': {
        'path': '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/lightdock_cluster_metrics.csv',
        'cluster_col': 'name',
        'label_col': 'category'
    },
    'Foldseek': {
        'path': '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/foldseek_cluster_metrics_enhanced.csv',
        'cluster_col': 'name',
        'label_col': 'category'
    },
    'Promoter': {
        'path': '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/promoter_similarity_results.csv',
        'cluster_col': 'group_name',
        'label_col': 'dataset_group'
    },
    'E2P2': {
        'path': '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/kegg_random/e2p2_feature_extraction_detailed.csv',
        'cluster_col': 'cluster_name',
        'label_col': 'classification_label'
    }
}

# Model configurations
MODELS = [
    {
        'name': 'CladePP',
        'script': '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/cladepp_phylo_profiling/cladepp_ml_classification_model.py',
        'summary_file': '/groups/itay_mayrose/alongonda/desktop/cladepp_ml_classification_summary.csv',
        'weights_file': '/groups/itay_mayrose/alongonda/desktop/cladepp_multi_feature_weights.csv',
        'features': [
            'mean_cladepp_score',
            'weighted_cladepp_score',
            'positive_correlation_ratio',
            'cladepp_multi_clade_high',
            'cladepp_multi_clade_medium',
            'cladepp_conservation_consistency',
            'cladepp_max_pair_score'
        ]
    },
    {
        'name': 'Docking',
        'script': '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/docking_matching/docking_ml_classification_model.py',
        'summary_file': '/groups/itay_mayrose/alongonda/desktop/docking_ml_classification_summary.csv',
        'weights_file': '/groups/itay_mayrose/alongonda/desktop/docking_multi_feature_weights.csv',
        'features': [
            'fraction_weak_binders',
            'q75_score',
            'fraction_strong_binders',
            'z_score',
            'max_score',
            'median_score_non_self',
            'mean_score_non_self'
        ]
    },
    {
        'name': 'Foldseek',
        'script': '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/foldseek_matching/foldseek_ml_classification_model.py',
        'summary_file': '/groups/itay_mayrose/alongonda/desktop/foldseek_ml_classification_summary.csv',
        'weights_file': '/groups/itay_mayrose/alongonda/desktop/foldseek_multi_feature_weights.csv',
        'features': [
            'mean_score_non_self',
            'enrichment_score',
            'z_score',
            'effect_size',
            'foldseek_match_coverage'
        ]
    },
    {
        'name': 'Promoter',
        'script': '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/extract_promotor/promoter_ml_classification_model.py',
        'summary_file': '/groups/itay_mayrose/alongonda/desktop/promoter_ml_classification_summary.csv',
        'weights_file': '/groups/itay_mayrose/alongonda/desktop/promoter_multi_feature_weights.csv',
        'features': [
            'mean_proximal_similarity',
            'mean_distal_similarity',
            'mean_proximal_correlation',
            'mean_distal_correlation',
            'num_tfbs_types_found'
        ]
    },
    {
        'name': 'E2P2',
        'script': '/groups/itay_mayrose/alongonda/desktop/python_scripts/features/structural_tailoring_classigication/e2p2_tagging/e2p2_ml_classification_model.py',
        'summary_file': '/groups/itay_mayrose/alongonda/desktop/e2p2_classification_models_summary.csv',
        'weights_file': '/groups/itay_mayrose/alongonda/desktop/random_kegg_e2p2_multi_feature_weights.csv',
        'features': [
            'num_distinct_enzyme_classes',
            'num_distinct_enzyme_subclasses',
            'num_distinct_enzyme_families',
            'num_distinct_enzyme_subfamilies',
            'total_ec_numbers'
        ]
    }
]

def run_model(model_config):
    """Run a single ML classification model"""
    print(f"\n{'='*100}")
    print(f"Running {model_config['name']} Model")
    print(f"{'='*100}")
    
    try:
        result = subprocess.run(
            [sys.executable, model_config['script']],
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )
        
        if result.returncode == 0:
            print(f"✅ {model_config['name']} model completed successfully")
            if result.stdout:
                # Print last 50 lines of output
                lines = result.stdout.strip().split('\n')
                print('\n'.join(lines[-50:]))
            return True
        else:
            print(f"❌ {model_config['name']} model failed with return code {result.returncode}")
            if result.stderr:
                print("Error output:")
                print(result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print(f"❌ {model_config['name']} model timed out after 10 minutes")
        return False
    except Exception as e:
        print(f"❌ Error running {model_config['name']} model: {e}")
        return False

def load_summary_results():
    """Load all summary results from CSV files"""
    all_results = []
    
    for model in MODELS:
        summary_file = Path(model['summary_file'])
        if not summary_file.exists():
            print(f"⚠️  Summary file not found: {summary_file}")
            continue
        
        try:
            df = pd.read_csv(summary_file)
            df['category'] = model['name']
            
            # Standardize column names: rename test_* columns to remove prefix
            # E2P2 uses test_roc_auc, test_pr_auc, etc. instead of roc_auc, pr_auc
            column_mapping = {}
            for col in df.columns:
                if col.startswith('test_'):
                    new_col = col.replace('test_', '')
                    column_mapping[col] = new_col
            
            if column_mapping:
                df = df.rename(columns=column_mapping)
                print(f"  Standardized columns: {list(column_mapping.values())}")
            
            all_results.append(df)
            print(f"✅ Loaded {model['name']}: {len(df)} results")
        except Exception as e:
            print(f"❌ Error loading {model['name']}: {e}")
            continue
    
    if not all_results:
        print("❌ No results loaded!")
        return None
    
    combined_df = pd.concat(all_results, ignore_index=True)
    return combined_df

def load_and_merge_all_data():
    """Load all feature data from all categories and merge on cluster names"""
    print("\n" + "="*100)
    print("LOADING AND MERGING ALL FEATURE DATA")
    print("="*100)
    
    merged_df = None
    
    for model in MODELS:
        category = model['name']
        data_config = DATA_FILES[category]
        
        try:
            df = pd.read_csv(data_config['path'])
            
            # Filter promoter data for successful analyses
            if category == 'Promoter' and 'status' in df.columns:
                df = df[df['status'] == 'success'].copy()
                print(f"  Filtered to {len(df)} successful analyses")
            
            # Create label
            if data_config['label_col'] == 'name':
                df['label'] = df[data_config['label_col']].apply(lambda x: 0 if 'RANDOM' in str(x) else 1)
            elif data_config['label_col'] == 'category':
                df['label'] = df[data_config['label_col']].apply(lambda x: 0 if x == 'KEGG_Random' or 'RANDOM' in str(x) else 1)
            elif data_config['label_col'] == 'dataset_group':
                df['label'] = df[data_config['label_col']].apply(lambda x: 0 if x == 'Random' else 1)
            elif data_config['label_col'] == 'classification_label':
                df['label'] = df[data_config['label_col']].astype(int)
            else:
                raise ValueError(f"Unknown label column: {data_config['label_col']}")
            
            print(f"  Labels: {df['label'].value_counts().to_dict()}")
            
            # Extract cluster_size
            if category == 'Docking' or category == 'Foldseek':
                # These have 'n' column for cluster size
                if 'n' in df.columns:
                    df['cluster_size'] = df['n']
                else:
                    print(f"⚠️  Warning: 'n' column not found in {category}, skipping cluster_size")
            elif category == 'Promoter':
                # Promoter has 'num_promoters' column
                if 'num_promoters' in df.columns:
                    df['cluster_size'] = df['num_promoters']
                else:
                    print(f"⚠️  Warning: 'num_promoters' column not found in {category}, skipping cluster_size")
            elif category == 'CladePP':
                # CladePP might have 'total_tips' or we need to compute from name
                if 'total_tips' in df.columns:
                    df['cluster_size'] = df['total_tips']
                else:
                    print(f"⚠️  Warning: 'total_tips' column not found in {category}, skipping cluster_size")
            elif category == 'E2P2':
                # E2P2 doesn't have cluster size directly, might need to compute
                # For now, skip if not available
                if 'cluster_size' not in df.columns:
                    print(f"⚠️  Warning: cluster_size not available in {category}, will try to get from merged data")
            
            # Select columns
            cluster_col = data_config['cluster_col']
            cols_to_keep = [cluster_col, 'label']
            
            # Add cluster_size if available
            if 'cluster_size' in df.columns:
                cols_to_keep.append('cluster_size')
            
            # Add features
            for feat in model['features']:
                if feat in df.columns:
                    cols_to_keep.append(feat)
                else:
                    print(f"⚠️  Warning: Feature '{feat}' not found in {category} data")
            
            df = df[cols_to_keep]
            
            # Rename cluster column to standard name
            df = df.rename(columns={cluster_col: 'cluster_name'})
            
            # Add prefix to feature names to avoid collisions
            prefix = category.lower() + '_'
            feature_rename = {}
            for feat in model['features']:
                if feat in df.columns:
                    feature_rename[feat] = prefix + feat
            
            # Rename cluster_size only once (use first occurrence)
            if 'cluster_size' in df.columns and 'cluster_size' not in feature_rename:
                # Keep cluster_size without prefix (it's a general feature)
                pass
            
            df = df.rename(columns=feature_rename)
            
            print(f"\n{category:12s}: {len(df):6d} clusters, {len([c for c in df.columns if c.startswith(prefix)])} features")
            if 'cluster_size' in df.columns:
                print(f"             cluster_size available")
            
            # Merge
            if merged_df is None:
                merged_df = df
            else:
                # Merge on cluster_name only first, then verify labels match
                merged_df = merged_df.merge(df, on=['cluster_name'], how='inner', suffixes=('', '_new'))
                
                # Handle label columns - they should match, but if not, use the first one
                if 'label_new' in merged_df.columns:
                    # Check for mismatches
                    mismatches = (merged_df['label'] != merged_df['label_new']).sum()
                    if mismatches > 0:
                        print(f"  ⚠️  Warning: {mismatches} label mismatches found, using first label")
                    merged_df = merged_df.drop(columns=['label_new'])
                
                # If cluster_size exists in both, keep the first one (or fill from second)
                if 'cluster_size' in merged_df.columns:
                    # Already have cluster_size, keep it
                    pass
                elif 'cluster_size_new' in merged_df.columns:
                    merged_df['cluster_size'] = merged_df['cluster_size_new']
                    merged_df = merged_df.drop(columns=['cluster_size_new'])
        
        except Exception as e:
            print(f"❌ Error loading {category}: {e}")
            continue
    
    if merged_df is None:
        print("❌ No data loaded!")
        return None
    
    print(f"\n✅ Merged dataset: {len(merged_df)} clusters")
    print(f"   MGC:    {merged_df['label'].sum()}")
    print(f"   Random: {(merged_df['label']==0).sum()}")
    print(f"   Total features: {len([c for c in merged_df.columns if c not in ['cluster_name', 'label', 'cluster_size']])}")
    if 'cluster_size' in merged_df.columns:
        print(f"   cluster_size: available")
    
    return merged_df

def train_combined_model(df):
    """Train a single combined model using ALL features from all categories + cluster_size
    with hyperparameter tuning, detailed progress, and advanced training features"""
    print("\n" + "="*100)
    print("COMBINED MODEL - ALL FEATURES FROM ALL CATEGORIES")
    print("Advanced Training with Hyperparameter Tuning")
    print("="*100)
    
    # Collect all features
    all_features = []
    for model in MODELS:
        prefix = model['name'].lower() + '_'
        for feat in model['features']:
            prefixed_feat = prefix + feat
            if prefixed_feat in df.columns:
                all_features.append(prefixed_feat)
    
    # Add cluster_size if available
    if 'cluster_size' in df.columns:
        all_features.insert(0, 'cluster_size')  # Add at the beginning
    
    print(f"\n📊 FEATURE COLLECTION:")
    print(f"  Total features: {len(all_features)}")
    if 'cluster_size' in all_features:
        print(f"  ✓ cluster_size (general feature)")
    for model in MODELS:
        prefix = model['name'].lower() + '_'
        model_features = [f for f in all_features if f.startswith(prefix)]
        if model_features:
            print(f"  ✓ {model['name']}: {len(model_features)} features")
    
    # Initialize Ray cluster early
    if not ray.is_initialized():
        print(f"\n  🚀 Initializing Ray cluster...")
        ray.init(ignore_reinit_error=True, num_cpus=None)
        print(f"    ✓ Ray initialized")
        available_cpus = int(ray.available_resources().get('CPU', 1))
        print(f"    - Available CPUs: {available_cpus}")
        print(f"    - Available GPUs: {ray.available_resources().get('GPU', 0)}")
    else:
        available_cpus = int(ray.available_resources().get('CPU', 1))
    
    # Prepare data
    print(f"\n{'='*100}")
    print("STEP 1: DATA PREPARATION WITH RAY DATA")
    print(f"{'='*100}")
    
    # Prepare DataFrame with all features and labels
    data = df[['cluster_name'] + all_features + ['label']].dropna()

    # Separate BGC clusters (validation-only) based on cluster_name containing 'BGC'
    bgc_mask = data['cluster_name'].str.contains('BGC', case=False, na=False)
    bgc_data = data[bgc_mask].copy()
    train_data = data[~bgc_mask].copy()
    
    print(f"  ✓ Loaded {len(data)} total samples")
    print(f"    - Training pool (non-BGC): {len(train_data)} samples")
    print(f"        * MGC (positive): {train_data['label'].sum()} ({100*train_data['label'].sum()/len(train_data):.1f}%)")
    print(f"        * Random (negative): {(train_data['label']==0).sum()} ({100*(train_data['label']==0).sum()/len(train_data):.1f}%)")
    print(f"    - Validation set (BGC only): {len(bgc_data)} samples")
    print(f"        * MGC (positive): {bgc_data['label'].sum()} ({100*bgc_data['label'].sum()/max(len(bgc_data),1):.1f}%)")
    print(f"    - Feature matrix (training) shape: ({len(train_data)}, {len(all_features)})")
    print(f"    - Missing values (training): {train_data[all_features].isna().sum().sum()} (should be 0)")
    
    # Convert to Ray Dataset immediately
    print(f"\n  📦 Converting training data to Ray Dataset...")
    dataset = from_pandas(train_data)
    print(f"    ✓ Ray Dataset created: {dataset.count()} samples (training only, non-BGC)")
    
    # Stratified train/test split using Ray Data
    print(f"\n  Splitting data (70% train, 30% test) with stratification...")
    def stratified_split(dataset, test_size=0.3, seed=42):
        """Stratified split by grouping by label and splitting each group"""
        # Convert to pandas temporarily for filtering (Ray Data filter syntax varies)
        df = dataset.to_pandas()
        
        # Split by label groups
        label_0_df = df[df['label'] == 0]
        label_1_df = df[df['label'] == 1]
        
        # Convert back to Ray Datasets
        label_0 = from_pandas(label_0_df)
        label_1 = from_pandas(label_1_df)
        
        # Split each group
        train_0, test_0 = label_0.train_test_split(test_size=test_size, seed=seed)
        train_1, test_1 = label_1.train_test_split(test_size=test_size, seed=seed)
        
        # Combine
        train_dataset = train_0.union(train_1)
        test_dataset = test_0.union(test_1)
        
        return train_dataset, test_dataset
    
    train_dataset, test_dataset = stratified_split(dataset, test_size=0.3, seed=42)
    
    # Count samples in each split
    train_count = train_dataset.count()
    test_count = test_dataset.count()
    
    # Get label counts for reporting
    train_df_temp = train_dataset.to_pandas()
    test_df_temp = test_dataset.to_pandas()
    train_mgc = train_df_temp['label'].sum()
    train_random = (train_df_temp['label'] == 0).sum()
    test_mgc = test_df_temp['label'].sum()
    test_random = (test_df_temp['label'] == 0).sum()
    
    print(f"    ✓ Training set: {train_count} samples ({train_mgc} MGC, {train_random} Random)")
    print(f"    ✓ Test set: {test_count} samples ({test_mgc} MGC, {test_random} Random)")
    
    # Distributed feature scaling with Ray Data
    print(f"\n  Scaling features using Ray Data (distributed StandardScaler)...")
    
    # Compute statistics from training data
    def compute_stats(batch):
        """Compute mean and std for each feature"""
        import pandas as pd
        feature_cols = [col for col in batch.columns if col in all_features]
        stats = {
            'mean': batch[feature_cols].mean().to_dict(),
            'std': batch[feature_cols].std().to_dict(),
            'count': len(batch)
        }
        return pd.DataFrame([stats])
    
    # Aggregate statistics from all batches
    train_stats_list = []
    for batch in train_dataset.iter_batches(batch_size=1000, batch_format="pandas"):
        stats = compute_stats(batch)
        train_stats_list.append(stats)
    
    # Combine statistics (weighted by count)
    total_count = sum(s['count'].iloc[0] for s in train_stats_list)
    mean_dict = {}
    std_dict = {}
    
    for feat in all_features:
        weighted_mean = sum(s['mean'].iloc[0].get(feat, 0) * s['count'].iloc[0] for s in train_stats_list) / total_count
        # For std, we need to compute pooled variance
        weighted_var = sum((s['std'].iloc[0].get(feat, 0)**2 + (s['mean'].iloc[0].get(feat, 0) - weighted_mean)**2) * s['count'].iloc[0] 
                          for s in train_stats_list) / total_count
        mean_dict[feat] = weighted_mean
        std_dict[feat] = np.sqrt(weighted_var) if weighted_var > 0 else 1.0
    
    print(f"    ✓ Statistics computed: mean/std for {len(all_features)} features")
    
    # Apply scaling using map_batches
    def scale_batch(batch, mean_dict, std_dict):
        """Scale features in a batch"""
        import pandas as pd
        batch = batch.copy()
        for feat in all_features:
            if feat in batch.columns:
                mean_val = mean_dict.get(feat, 0)
                std_val = std_dict.get(feat, 1.0)
                if std_val > 0:
                    batch[feat] = (batch[feat] - mean_val) / std_val
                else:
                    batch[feat] = batch[feat] - mean_val
        return batch
    
    # Apply scaling using map_batches (Ray Data will automatically parallelize)
    # Don't specify num_cpus here - let Ray Data handle parallelism automatically
    # This allows multiple trials to run in parallel
    train_dataset_scaled = train_dataset.map_batches(
        scale_batch, 
        fn_kwargs={'mean_dict': mean_dict, 'std_dict': std_dict},
        batch_format="pandas"
    )
    test_dataset_scaled = test_dataset.map_batches(
        scale_batch,
        fn_kwargs={'mean_dict': mean_dict, 'std_dict': std_dict},
        batch_format="pandas"
    )
    
    print(f"    ✓ Scaling applied to train and test datasets")
    
    # Verify scaling (convert to pandas for verification)
    train_scaled_df = train_dataset_scaled.to_pandas()
    test_scaled_df = test_dataset_scaled.to_pandas()
    train_scaled_features = train_scaled_df[all_features].values
    test_scaled_features = test_scaled_df[all_features].values
    
    print(f"    ✓ Training set scaled: mean={train_scaled_features.mean():.6f}, std={train_scaled_features.std():.6f}")
    print(f"    ✓ Test set scaled: mean={test_scaled_features.mean():.6f}, std={test_scaled_features.std():.6f}")
    
    # Hyperparameter tuning with Ray Tune
    print(f"\n{'='*100}")
    print("STEP 2: HYPERPARAMETER TUNING WITH RAY TUNE")
    print(f"{'='*100}")
    
    # Use Ray Tune for advanced hyperparameter tuning
    from sklearn.model_selection import cross_val_score
    
    print(f"  Using Ray Tune with Optuna search algorithm")
    print(f"  Features: Bayesian optimization, early stopping, distributed tuning")
    
    # Define search space for Ray Tune
    search_space = {
        "C": tune.loguniform(1e-3, 1e3),  # Log-uniform distribution
        "solver": tune.choice(["lbfgs", "liblinear", "saga"]),
        "class_weight": tune.choice([None, "balanced", {0: 1, 1: 2}, {0: 1, 1: 3}]),
        "penalty": tune.choice(["l1", "l2"])
    }
    
    print(f"\n  Search space:")
    print(f"    - C: log-uniform [0.001, 1000]")
    print(f"    - solver: ['lbfgs', 'liblinear', 'saga']")
    print(f"    - class_weight: [None, 'balanced', custom weights]")
    print(f"    - penalty: ['l1', 'l2']")
    print(f"  Scoring metric: PR AUC (average_precision)")
    print(f"  Cross-validation: 5-fold stratified")
    print(f"  Search algorithm: Optuna (Bayesian optimization)")
    print(f"  Early stopping: ASHA scheduler")
    
    start_time = time.time()
    
    # Define objective function for Ray Tune
    def train_model_tune(config):
        """Objective function for Ray Tune - trains scikit-learn model with given config"""
        from sklearn.linear_model import LogisticRegression
        from sklearn.model_selection import cross_val_score
        
        # Get data from Ray Dataset (convert to pandas for sklearn)
        train_data_local = train_dataset_scaled.to_pandas()
        X_train_ray = train_data_local[all_features].values
        y_train_ray = train_data_local['label'].values
        
        # Handle solver-penalty compatibility
        solver = config["solver"]
        penalty = config["penalty"]
        
        # LBFGS only supports L2
        if solver == "lbfgs" and penalty == "l1":
            penalty = "l2"
        
        # Create model with config
        model = LogisticRegression(
            C=config["C"],
            solver=solver,
            penalty=penalty,
            class_weight=config["class_weight"],
            random_state=42,
            max_iter=2000,
            n_jobs=1  # Ray handles parallelism
        )
        
        # 5-fold cross-validation
        cv_scores = cross_val_score(
            model, X_train_ray, y_train_ray,
            cv=5, scoring='average_precision', n_jobs=1
        )
        
        # Report metrics to Ray Tune
        tune.report({
            "mean_pr_auc": float(cv_scores.mean()),
            "std_pr_auc": float(cv_scores.std()),
            "min_pr_auc": float(cv_scores.min()),
            "max_pr_auc": float(cv_scores.max())
        })
    
    # Initialize Optuna search
    optuna_search = OptunaSearch(metric="mean_pr_auc", mode="max")
    
    # ASHA scheduler for early stopping
    scheduler = ASHAScheduler(
        metric="mean_pr_auc",
        mode="max",
        max_t=100,  # Maximum iterations
        grace_period=10,  # Minimum iterations before stopping
        reduction_factor=2
    )
    
    # Get available CPUs for parallel trials (if not already set)
    if 'available_cpus' not in locals():
        available_cpus = int(ray.available_resources().get('CPU', 1))
    # Use most CPUs for parallel trials (leave a few for system)
    num_parallel_trials = max(1, available_cpus - 2)
    
    print(f"\n  🔍 Running Ray Tune hyperparameter search...")
    print(f"    - Number of trials: 100 (with early stopping)")
    print(f"    - Parallel trials: {num_parallel_trials} (using {available_cpus} available CPUs)")
    
    # Create Tuner
    tuner = tune.Tuner(
        train_model_tune,
        tune_config=tune.TuneConfig(
            search_alg=optuna_search,
            scheduler=scheduler,
            num_samples=100,  # Number of trials
            # metric and mode are specified in scheduler, don't duplicate here
        ),
        param_space=search_space,
        run_config=tune.RunConfig(
            name="combined_model_tuning",
            stop={"training_iteration": 100},
            verbose=1,
        ),
    )
    
    # Run tuning
    results = tuner.fit()
    
    tuning_time = time.time() - start_time
    print(f"\n  ✓ Hyperparameter tuning completed in {tuning_time:.1f} seconds")
    
    # Get best result
    best_result = results.get_best_result(metric="mean_pr_auc", mode="max")
    best_params = best_result.config
    best_cv_score = best_result.metrics["mean_pr_auc"]
    
    print(f"\n  🏆 BEST PARAMETERS (from {len(results)} trials):")
    for param, value in best_params.items():
        print(f"    - {param}: {value}")
    print(f"  Best CV PR AUC: {best_cv_score:.6f} ± {best_result.metrics.get('std_pr_auc', 0):.6f}")
    print(f"  Min PR AUC: {best_result.metrics.get('min_pr_auc', 0):.6f}")
    print(f"  Max PR AUC: {best_result.metrics.get('max_pr_auc', 0):.6f}")
    
    # Create best model with best parameters
    solver = best_params["solver"]
    penalty = best_params["penalty"]
    if solver == "lbfgs" and penalty == "l1":
        penalty = "l2"
    
    best_model = LogisticRegression(
        C=best_params["C"],
        solver=solver,
        penalty=penalty,
        class_weight=best_params["class_weight"],
        random_state=42,
        max_iter=2000,
        n_jobs=-1
    )
    
    # Save tuning results
    tuning_results_df = pd.DataFrame([
        {**r.config, **r.metrics} for r in results
    ])
    combined_output_dir = "/groups/itay_mayrose/alongonda/desktop/mibig_validate/combined"
    os.makedirs(combined_output_dir, exist_ok=True)
    tuning_file = os.path.join(combined_output_dir, "combined_model_ray_tune_results.csv")
    tuning_results_df.to_csv(tuning_file, index=False)
    print(f"  💾 Ray Tune results saved to: {tuning_file}")
    
    # Training with best model
    print(f"\n{'='*100}")
    print("STEP 3: TRAINING BEST MODEL")
    print(f"{'='*100}")
    
    # Convert Ray Dataset to pandas/numpy for sklearn training
    train_scaled_df = train_dataset_scaled.to_pandas()
    test_scaled_df = test_dataset_scaled.to_pandas()
    X_train_scaled = train_scaled_df[all_features].values
    y_train = train_scaled_df['label'].values
    X_test_scaled = test_scaled_df[all_features].values
    y_test = test_scaled_df['label'].values
    
    print(f"  Training final model with best parameters...")
    train_start = time.time()
    best_model.fit(X_train_scaled, y_train)
    train_time = time.time() - train_start
    print(f"    ✓ Model trained in {train_time:.2f} seconds")
    print(f"    ✓ Converged: {best_model.n_iter_[0]} iterations")
    
    # Get weights
    weights = best_model.coef_[0]
    intercept = best_model.intercept_[0]
    print(f"    ✓ Model coefficients: {len(weights)} features")
    print(f"    ✓ Intercept: {intercept:.6f}")
    print(f"    ✓ Non-zero coefficients: {(np.abs(weights) > 1e-6).sum()}")
    
    # Learning curves
    print(f"\n  📈 Computing learning curves...")
    train_sizes, train_scores, val_scores = learning_curve(
        best_model, X_train_scaled, y_train, cv=5, 
        scoring='average_precision', n_jobs=-1,
        train_sizes=np.linspace(0.1, 1.0, 10), verbose=0
    )
    print(f"    ✓ Learning curves computed")
    print(f"    Final training score: {train_scores[-1].mean():.6f} ± {train_scores[-1].std():.6f}")
    print(f"    Final validation score: {val_scores[-1].mean():.6f} ± {val_scores[-1].std():.6f}")
    
    # Validation curve for C parameter
    print(f"\n  📊 Computing validation curve for C parameter...")
    C_range = np.logspace(-3, 3, 7)
    val_curve_params = {
        'solver': best_params['solver'],
        'class_weight': best_params['class_weight'],
        'random_state': 42,
        'max_iter': 2000
    }
    # Only add penalty if solver supports it
    if best_params['solver'] in ['liblinear', 'saga']:
        val_curve_params['penalty'] = best_params['penalty']
    
    train_scores_cv, val_scores_cv = validation_curve(
        LogisticRegression(**val_curve_params),
        X_train_scaled, y_train, param_name='C', param_range=C_range,
        cv=5, scoring='average_precision', n_jobs=-1, verbose=0
    )
    print(f"    ✓ Validation curve computed")
    best_c_idx = np.argmax(val_scores_cv.mean(axis=1))
    print(f"    Best C from validation curve: {C_range[best_c_idx]:.3f} (PR AUC: {val_scores_cv.mean(axis=1)[best_c_idx]:.6f})")
    
    # Predictions
    print(f"\n{'='*100}")
    print("STEP 4: MODEL EVALUATION")
    print(f"{'='*100}")
    
    print(f"  Making predictions on test set (non-BGC)...")
    y_pred_proba = best_model.predict_proba(X_test_scaled)[:, 1]
    print(f"    ✓ Predictions completed")
    print(f"    Probability range (test): [{y_pred_proba.min():.4f}, {y_pred_proba.max():.4f}]")
    print(f"    Mean probability (test): {y_pred_proba.mean():.4f}")
    
    # ROC curve
    fpr, tpr, roc_thresholds = roc_curve(y_test, y_pred_proba)
    roc_auc = auc(fpr, tpr)
    
    # Precision-Recall curve
    precision, recall, pr_thresholds = precision_recall_curve(y_test, y_pred_proba)
    pr_auc = auc(recall, precision)
    
    # Find F1-optimal threshold
    f1_scores = 2 * (precision[:-1] * recall[:-1]) / (precision[:-1] + recall[:-1] + 1e-10)
    f1_optimal_idx = np.argmax(f1_scores)
    f1_optimal_threshold = pr_thresholds[f1_optimal_idx]
    
    print(f"\n  📊 Performance Metrics:")
    print(f"    ROC AUC: {roc_auc:.6f}")
    print(f"    PR AUC:  {pr_auc:.6f}")
    print(f"    F1-optimal threshold: {f1_optimal_threshold:.4f}")
    
    # Evaluate at F1-optimal threshold
    y_pred = (y_pred_proba >= f1_optimal_threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
    
    precision_val = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall_val = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = f1_score(y_test, y_pred)
    mcc = matthews_corrcoef(y_test, y_pred)
    
    print(f"\n  🎯 Performance at F1-optimal threshold (TEST SET, non-BGC):")
    print(f"    Precision: {precision_val:.6f}")
    print(f"    Recall:    {recall_val:.6f}")
    print(f"    F1 Score:  {f1:.6f}")
    print(f"    MCC:       {mcc:.6f}")
    print(f"\n  📋 Confusion Matrix:")
    print(f"    True Positives (TP):  {tp:4d}")
    print(f"    False Positives (FP): {fp:4d}")
    print(f"    True Negatives (TN):  {tn:5d}")
    print(f"    False Negatives (FN): {fn:4d}")

    # Evaluate on validation set (BGC-only) using same threshold
    print(f"\n  ✅ Evaluating on BGC VALIDATION SET (held-out, all BGC clusters)...")
    if len(bgc_data) > 0:
        # Scale BGC features using training statistics
        bgc_scaled = bgc_data.copy()
        for feat in all_features:
            if feat in bgc_scaled.columns:
                mean_val = mean_dict.get(feat, 0)
                std_val = std_dict.get(feat, 1.0)
                if std_val > 0:
                    bgc_scaled[feat] = (bgc_scaled[feat] - mean_val) / std_val
                else:
                    bgc_scaled[feat] = bgc_scaled[feat] - mean_val

        X_val_scaled = bgc_scaled[all_features].values
        y_val = bgc_scaled['label'].values

        y_pred_proba_val = best_model.predict_proba(X_val_scaled)[:, 1]
        print(f"    ✓ Predictions completed on BGC validation set")
        print(f"    Probability range (BGC): [{y_pred_proba_val.min():.4f}, {y_pred_proba_val.max():.4f}]")
        print(f"    Mean probability (BGC): {y_pred_proba_val.mean():.4f}")

        # ROC & PR on validation set
        try:
            fpr_val, tpr_val, _ = roc_curve(y_val, y_pred_proba_val)
            roc_auc_val = auc(fpr_val, tpr_val)
        except ValueError:
            roc_auc_val = np.nan

        precision_val_curve, recall_val_curve, _ = precision_recall_curve(y_val, y_pred_proba_val)
        pr_auc_val = auc(recall_val_curve, precision_val_curve)

        # Evaluate validation set at F1-optimal threshold from test set
        y_pred_val = (y_pred_proba_val >= f1_optimal_threshold).astype(int)
        tn_val, fp_val, fn_val, tp_val = confusion_matrix(y_val, y_pred_val).ravel()

        val_precision = tp_val / (tp_val + fp_val) if (tp_val + fp_val) > 0 else 0
        val_recall = tp_val / (tp_val + fn_val) if (tp_val + fn_val) > 0 else 0
        f1_val = f1_score(y_val, y_pred_val)
        mcc_val = matthews_corrcoef(y_val, y_pred_val) if (tp_val + fp_val + tn_val + fn_val) > 0 else 0

        print(f"\n  📊 BGC VALIDATION PERFORMANCE (at test-set F1-optimal threshold):")
        print(f"    BGC count: {len(y_val)}")
        print(f"    ROC AUC: {roc_auc_val:.6f}" if not np.isnan(roc_auc_val) else "    ROC AUC: NaN (only positive class present)")
        print(f"    PR AUC:  {pr_auc_val:.6f}")
        print(f"    Precision: {val_precision:.6f}")
        print(f"    Recall:    {val_recall:.6f}  ({tp_val}/{len(y_val)} BGC correctly detected, {fn_val}/{len(y_val)} missed)")
        print(f"    F1 Score:  {f1_val:.6f}")
        print(f"    MCC:       {mcc_val:.6f}")
        print(f"\n  📋 BGC Confusion Matrix:")
        print(f"    True Positives (TP):  {tp_val:4d}")
        print(f"    False Positives (FP): {fp_val:4d}")
        print(f"    True Negatives (TN):  {tn_val:5d}")
        print(f"    False Negatives (FN): {fn_val:4d}")
    else:
        roc_auc_val = np.nan
        pr_auc_val = np.nan
        val_precision = np.nan
        val_recall = np.nan
        f1_val = np.nan
        mcc_val = np.nan
        tp_val = 0
        fn_val = 0
    
    # Cross-validation on full dataset (combine train and test)
    print(f"\n  🔄 Running 5-fold cross-validation on full dataset...")
    cv_start = time.time()
    # Combine train and test datasets for full cross-validation
    full_dataset_scaled = train_dataset_scaled.union(test_dataset_scaled)
    full_df = full_dataset_scaled.to_pandas()
    X_full = full_df[all_features].values
    y_full = full_df['label'].values
    cv_scores = cross_validate(best_model, X_full, y_full, cv=5, 
                               scoring=['roc_auc', 'average_precision', 'f1', 'precision', 'recall'], 
                               return_train_score=True, n_jobs=-1)
    cv_time = time.time() - cv_start
    print(f"    ✓ Cross-validation completed in {cv_time:.1f} seconds")
    
    print(f"\n  📊 Cross-Validation Results:")
    print(f"    ROC AUC: {cv_scores['test_roc_auc'].mean():.6f} ± {cv_scores['test_roc_auc'].std():.6f}")
    print(f"    PR AUC:  {cv_scores['test_average_precision'].mean():.6f} ± {cv_scores['test_average_precision'].std():.6f}")
    print(f"    F1:      {cv_scores['test_f1'].mean():.6f} ± {cv_scores['test_f1'].std():.6f}")
    print(f"    Precision: {cv_scores['test_precision'].mean():.6f} ± {cv_scores['test_precision'].std():.6f}")
    print(f"    Recall:    {cv_scores['test_recall'].mean():.6f} ± {cv_scores['test_recall'].std():.6f}")
    
    # Feature importance
    print(f"\n{'='*100}")
    print("STEP 5: FEATURE IMPORTANCE ANALYSIS")
    print(f"{'='*100}")
    
    feature_importance = sorted(zip(all_features, np.abs(weights)), key=lambda x: x[1], reverse=True)
    print(f"\n  🏆 ALL FEATURES BY IMPORTANCE (Ranked):")
    print("-" * 100)
    print(f"{'Rank':<6} {'Feature':<60} {'Weight':<12} {'Importance':<12}")
    print("-" * 100)
    for rank, (feat, imp) in enumerate(feature_importance, 1):
        weight = weights[all_features.index(feat)]
        sign = '+' if weight >= 0 else '-'
        print(f"{rank:<6} {feat:<60} {sign}{abs(weight):<11.6f} {imp:<12.6f}")
    
    # Save weights
    weights_df = pd.DataFrame({
        'feature': all_features,
        'weight': weights,
        'abs_importance': np.abs(weights)
    }).sort_values('abs_importance', ascending=False)
    
    combined_output_dir = "/groups/itay_mayrose/alongonda/desktop/mibig_validate/combined"
    os.makedirs(combined_output_dir, exist_ok=True)
    weights_file = os.path.join(combined_output_dir, "combined_all_features_weights.csv")
    weights_df.to_csv(weights_file, index=False)
    print(f"\n  💾 Combined model weights saved to: {weights_file}")
    
    # Save hyperparameter tuning results (if not already saved by Ray Tune)
    # Ray Tune results are saved above
    
    total_time = time.time() - start_time
    print(f"\n{'='*100}")
    print(f"✅ TRAINING COMPLETE - Total time: {total_time:.1f} seconds")
    print(f"{'='*100}")
    
    return {
        'roc_auc': roc_auc,
        'pr_auc': pr_auc,
        'precision': precision_val,
        'recall': recall_val,
        'val_roc_auc': roc_auc_val,
        'val_pr_auc': pr_auc_val,
        'val_precision': val_precision,
        'val_recall': val_recall,
        'val_f1': f1_val,
        'val_mcc': mcc_val,
        'val_tp': tp_val,
        'val_fn': fn_val,
        'val_size': int(len(bgc_data)),
        'f1': f1,
        'mcc': mcc,
        'best_params': best_params,
        'best_cv_score': best_cv_score,
        'cv_roc_auc': cv_scores['test_roc_auc'].mean(),
        'cv_roc_auc_std': cv_scores['test_roc_auc'].std(),
        'cv_pr_auc': cv_scores['test_average_precision'].mean(),
        'cv_pr_auc_std': cv_scores['test_average_precision'].std(),
        'cv_f1': cv_scores['test_f1'].mean(),
        'cv_f1_std': cv_scores['test_f1'].std(),
        'cv_precision': cv_scores['test_precision'].mean(),
        'cv_precision_std': cv_scores['test_precision'].std(),
        'cv_recall': cv_scores['test_recall'].mean(),
        'cv_recall_std': cv_scores['test_recall'].std(),
        'n_features': len(all_features),
        'n_samples': len(data),
        'training_time': total_time
    }

def create_comprehensive_analysis():
    """Create comprehensive analysis of all ML models"""
    print("\n" + "="*100)
    print("COMPREHENSIVE ML CLASSIFICATION ANALYSIS")
    print("All 5 Feature Categories")
    print("="*100)
    
    # Step 1: Run all models
    print("\n" + "="*100)
    print("STEP 1: RUNNING ALL ML MODELS")
    print("="*100)
    
    results = {}
    for model in MODELS:
        success = run_model(model)
        results[model['name']] = success
    
    # Step 2: Load all results
    print("\n" + "="*100)
    print("STEP 2: LOADING RESULTS")
    print("="*100)
    
    combined_df = load_summary_results()
    if combined_df is None:
        print("❌ Cannot proceed without results")
        return None
    
    # Step 3: Individual Feature Analysis
    print("\n" + "="*100)
    print("STEP 3: INDIVIDUAL FEATURE PERFORMANCE")
    print("="*100)
    
    individual_df = combined_df[combined_df['model_type'] == 'individual'].copy()
    
    print(f"\nTotal individual features analyzed: {len(individual_df)}")
    print(f"Categories: {individual_df['category'].unique()}")
    
    # Rank all features by PR AUC (handle NaN values)
    if 'pr_auc' in individual_df.columns:
        individual_df_sorted = individual_df.sort_values('pr_auc', ascending=False, na_position='last')
    else:
        individual_df_sorted = individual_df.sort_values('f1', ascending=False, na_position='last')
    
    print(f"\n{'='*120}")
    print("TOP 20 INDIVIDUAL FEATURES (Ranked by PR AUC)")
    print(f"{'='*120}")
    print(f"{'Rank':<6} {'Category':<12} {'Feature':<45} {'ROC AUC':<10} {'PR AUC':<10} {'Precision':<10} {'Recall':<10} {'F1':<10}")
    print("-" * 120)
    
    for idx, (_, row) in enumerate(individual_df_sorted.head(20).iterrows(), 1):
        roc_auc_str = f"{row['roc_auc']:.4f}" if pd.notna(row.get('roc_auc')) else "N/A"
        pr_auc_str = f"{row['pr_auc']:.4f}" if pd.notna(row.get('pr_auc')) else "N/A"
        precision_str = f"{row['precision']:.4f}" if pd.notna(row.get('precision')) else "N/A"
        recall_str = f"{row['recall']:.4f}" if pd.notna(row.get('recall')) else "N/A"
        f1_str = f"{row['f1']:.4f}" if pd.notna(row.get('f1')) else "N/A"
        print(f"{idx:<6} {row['category']:<12} {str(row['feature']):<45} {roc_auc_str:<10} {pr_auc_str:<10} {precision_str:<10} {recall_str:<10} {f1_str:<10}")
    
    # Step 4: Category-wise Analysis
    print(f"\n{'='*120}")
    print("PERFORMANCE BY CATEGORY")
    print(f"{'='*120}")
    
    for category in sorted(individual_df['category'].unique()):
        cat_df = individual_df[individual_df['category'] == category]
        if 'pr_auc' in cat_df.columns:
            cat_df_sorted = cat_df.sort_values('pr_auc', ascending=False, na_position='last')
        else:
            cat_df_sorted = cat_df.sort_values('f1', ascending=False, na_position='last')
        
        print(f"\n📊 {category.upper()} ({len(cat_df)} features)")
        print("-" * 120)
        print(f"{'Feature':<50} {'ROC AUC':<10} {'PR AUC':<10} {'Precision':<10} {'Recall':<10} {'F1':<10}")
        print("-" * 120)
        
        for _, row in cat_df_sorted.iterrows():
            roc_auc_str = f"{row['roc_auc']:.4f}" if pd.notna(row.get('roc_auc')) else "N/A"
            pr_auc_str = f"{row['pr_auc']:.4f}" if pd.notna(row.get('pr_auc')) else "N/A"
            precision_str = f"{row['precision']:.4f}" if pd.notna(row.get('precision')) else "N/A"
            recall_str = f"{row['recall']:.4f}" if pd.notna(row.get('recall')) else "N/A"
            f1_str = f"{row['f1']:.4f}" if pd.notna(row.get('f1')) else "N/A"
            print(f"{str(row['feature']):<50} {roc_auc_str:<10} {pr_auc_str:<10} {precision_str:<10} {recall_str:<10} {f1_str:<10}")
        
        # Category statistics
        print(f"\n  Category Statistics:")
        if 'pr_auc' in cat_df.columns and cat_df['pr_auc'].notna().any():
            best_pr_idx = cat_df['pr_auc'].idxmax()
            print(f"    Best PR AUC: {cat_df['pr_auc'].max():.4f} ({cat_df.loc[best_pr_idx, 'feature']})")
            print(f"    Mean PR AUC: {cat_df['pr_auc'].mean():.4f} ± {cat_df['pr_auc'].std():.4f}")
        if 'roc_auc' in cat_df.columns and cat_df['roc_auc'].notna().any():
            best_roc_idx = cat_df['roc_auc'].idxmax()
            print(f"    Best ROC AUC: {cat_df['roc_auc'].max():.4f} ({cat_df.loc[best_roc_idx, 'feature']})")
            print(f"    Mean ROC AUC: {cat_df['roc_auc'].mean():.4f} ± {cat_df['roc_auc'].std():.4f}")
    
    # Step 5: Multi-Feature Model Analysis
    print(f"\n{'='*120}")
    print("MULTI-FEATURE MODEL PERFORMANCE")
    print(f"{'='*120}")
    
    multi_df = combined_df[combined_df['model_type'] == 'multi_feature'].copy()
    if 'pr_auc' in multi_df.columns:
        multi_df_sorted = multi_df.sort_values('pr_auc', ascending=False, na_position='last')
    else:
        multi_df_sorted = multi_df.sort_values('f1', ascending=False, na_position='last')
    
    print(f"{'Category':<15} {'ROC AUC':<12} {'PR AUC':<12} {'Precision':<12} {'Recall':<12} {'F1':<12} {'MCC':<12}")
    print("-" * 120)
    
    for _, row in multi_df_sorted.iterrows():
        roc_auc_str = f"{row['roc_auc']:.4f}" if pd.notna(row.get('roc_auc')) else "N/A"
        pr_auc_str = f"{row['pr_auc']:.4f}" if pd.notna(row.get('pr_auc')) else "N/A"
        precision_str = f"{row['precision']:.4f}" if pd.notna(row.get('precision')) else "N/A"
        recall_str = f"{row['recall']:.4f}" if pd.notna(row.get('recall')) else "N/A"
        f1_str = f"{row['f1']:.4f}" if pd.notna(row.get('f1')) else "N/A"
        mcc_str = f"{row['mcc']:.4f}" if pd.notna(row.get('mcc')) else "N/A"
        print(f"{row['category']:<15} {roc_auc_str:<12} {pr_auc_str:<12} {precision_str:<12} {recall_str:<12} {f1_str:<12} {mcc_str:<12}")
    
    # Step 6: Feature Importance Analysis
    print(f"\n{'='*120}")
    print("FEATURE IMPORTANCE IN MULTI-FEATURE MODELS")
    print(f"{'='*120}")
    
    for model in MODELS:
        weights_file = Path(model['weights_file'])
        if not weights_file.exists():
            continue
        
        try:
            weights_df = pd.read_csv(weights_file)
            weights_df = weights_df.sort_values('abs_importance', ascending=False)
            
            print(f"\n📊 {model['name'].upper()} - Top Features by Importance")
            print("-" * 100)
            print(f"{'Feature':<50} {'Weight':<12} {'Abs Importance':<15}")
            print("-" * 100)
            
            for _, row in weights_df.head(10).iterrows():
                sign = '+' if row['weight'] >= 0 else '-'
                print(f"{str(row['feature']):<50} {sign}{abs(row['weight']):<11.6f} {row['abs_importance']:<15.6f}")
                
        except Exception as e:
            print(f"⚠️  Error loading weights for {model['name']}: {e}")
    
    # Step 7: Summary Statistics
    print(f"\n{'='*120}")
    print("OVERALL SUMMARY STATISTICS")
    print(f"{'='*120}")
    
    print(f"\nIndividual Features:")
    print(f"  Total features analyzed: {len(individual_df)}")
    if 'roc_auc' in individual_df.columns and individual_df['roc_auc'].notna().any():
        print(f"  Mean ROC AUC: {individual_df['roc_auc'].mean():.4f} ± {individual_df['roc_auc'].std():.4f}")
        print(f"  Best ROC AUC: {individual_df['roc_auc'].max():.4f}")
    if 'pr_auc' in individual_df.columns and individual_df['pr_auc'].notna().any():
        print(f"  Mean PR AUC: {individual_df['pr_auc'].mean():.4f} ± {individual_df['pr_auc'].std():.4f}")
        print(f"  Best PR AUC: {individual_df['pr_auc'].max():.4f}")
    if 'f1' in individual_df.columns and individual_df['f1'].notna().any():
        print(f"  Mean F1: {individual_df['f1'].mean():.4f} ± {individual_df['f1'].std():.4f}")
        print(f"  Best F1: {individual_df['f1'].max():.4f}")
    
    print(f"\nMulti-Feature Models:")
    print(f"  Total models: {len(multi_df)}")
    if 'roc_auc' in multi_df.columns and multi_df['roc_auc'].notna().any():
        print(f"  Mean ROC AUC: {multi_df['roc_auc'].mean():.4f} ± {multi_df['roc_auc'].std():.4f}")
        print(f"  Best ROC AUC: {multi_df['roc_auc'].max():.4f}")
    if 'pr_auc' in multi_df.columns and multi_df['pr_auc'].notna().any():
        print(f"  Mean PR AUC: {multi_df['pr_auc'].mean():.4f} ± {multi_df['pr_auc'].std():.4f}")
        print(f"  Best PR AUC: {multi_df['pr_auc'].max():.4f}")
    if 'f1' in multi_df.columns and multi_df['f1'].notna().any():
        print(f"  Mean F1: {multi_df['f1'].mean():.4f} ± {multi_df['f1'].std():.4f}")
        print(f"  Best F1: {multi_df['f1'].max():.4f}")
    
    # Step 8: Combined Model (All Features Together)
    print(f"\n{'='*120}")
    print("STEP 8: COMBINED MODEL - ALL FEATURES FROM ALL CATEGORIES")
    print(f"{'='*120}")
    
    merged_data = load_and_merge_all_data()
    if merged_data is not None:
        combined_result = train_combined_model(merged_data)
        
        print(f"\n{'='*120}")
        print("COMBINED MODEL SUMMARY")
        print(f"{'='*120}")
        print(f"Features: {combined_result['n_features']}")
        print(f"Samples: {combined_result['n_samples']}")
        print(f"Training time: {combined_result.get('training_time', 0):.1f} seconds")
        print(f"\nBest Hyperparameters:")
        for param, value in combined_result.get('best_params', {}).items():
            print(f"  - {param}: {value}")
        print(f"Best CV PR AUC: {combined_result.get('best_cv_score', 0):.6f}")
        print(f"\nPerformance Metrics:")
        print(f"ROC AUC: {combined_result['roc_auc']:.4f} (CV: {combined_result['cv_roc_auc']:.4f} ± {combined_result['cv_roc_auc_std']:.4f})")
        print(f"PR AUC:  {combined_result['pr_auc']:.4f} (CV: {combined_result['cv_pr_auc']:.4f} ± {combined_result['cv_pr_auc_std']:.4f})")
        print(f"Precision: {combined_result['precision']:.4f} (CV: {combined_result['cv_precision']:.4f} ± {combined_result['cv_precision_std']:.4f})")
        print(f"Recall:    {combined_result['recall']:.4f} (CV: {combined_result['cv_recall']:.4f} ± {combined_result['cv_recall_std']:.4f})")
        print(f"F1:        {combined_result['f1']:.4f} (CV: {combined_result['cv_f1']:.4f} ± {combined_result['cv_f1_std']:.4f})")
        print(f"MCC:       {combined_result['mcc']:.4f}")
    else:
        print("❌ Could not load merged data for combined model")
        combined_result = None
    
    # Step 9: Save comprehensive results
    output_file = "/groups/itay_mayrose/alongonda/desktop/comprehensive_ml_analysis_results.csv"
    combined_df.to_csv(output_file, index=False)
    print(f"\n💾 Comprehensive results saved to: {output_file}")
    
    # Save top features
    top_features_file = "/groups/itay_mayrose/alongonda/desktop/top_individual_features.csv"
    individual_df_sorted.to_csv(top_features_file, index=False)
    print(f"💾 Top individual features saved to: {top_features_file}")
    
    # Save multi-feature summary
    multi_summary_file = "/groups/itay_mayrose/alongonda/desktop/multi_feature_models_summary.csv"
    multi_df_sorted.to_csv(multi_summary_file, index=False)
    print(f"💾 Multi-feature models summary saved to: {multi_summary_file}")
    
    # Save combined model summary
    if combined_result:
        combined_output_dir = "/groups/itay_mayrose/alongonda/desktop/mibig_validate/combined"
        os.makedirs(combined_output_dir, exist_ok=True)
        combined_summary_file = os.path.join(combined_output_dir, "combined_all_features_summary.csv")
        combined_summary_df = pd.DataFrame([combined_result])
        combined_summary_df.to_csv(combined_summary_file, index=False)
        print(f"💾 Combined model summary saved to: {combined_summary_file}")
    
    return combined_df

if __name__ == "__main__":
    results_df = create_comprehensive_analysis()

