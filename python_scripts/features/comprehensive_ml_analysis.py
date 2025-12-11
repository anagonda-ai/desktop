#!/usr/bin/env python3
"""
Comprehensive ML Analysis - All 5 Feature Categories
Runs all ML classification models and creates comprehensive analysis
"""

import subprocess
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

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
            'similarity_score',
            'mean_proximal_similarity',
            'mean_distal_similarity',
            'correlation_score',
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
    
    # Step 8: Save comprehensive results
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
    
    return combined_df

if __name__ == "__main__":
    results_df = create_comprehensive_analysis()

