#!/usr/bin/env python3
"""
Docking Feature ML Classification Model - Individual & Multi-feature
"""

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (roc_curve, auc, precision_recall_curve, confusion_matrix, 
                            classification_report, matthews_corrcoef, f1_score)
from sklearn.model_selection import train_test_split, cross_validate
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

DOCKING_FEATURES = [
    'mean_score_non_self',
    'enrichment_score',
    'z_score',
    'effect_size'
]

FEATURE_DESCRIPTIONS = {
    'mean_score_non_self': 'Mean LightDock binding energy from all-vs-all protein-protein docking simulations, excluding self-interactions. More negative values indicate stronger predicted binding affinity. Aggregates pairwise binding potential across all protein pairs in the cluster',
    'enrichment_score': 'Fold-enrichment of observed mean binding energy vs. expected random baseline (computed from cluster size). Values with large absolute magnitude indicate non-random binding patterns. Formula: mean_score_non_self / expected_random',
    'z_score': 'Standardized score measuring how many standard deviations the cluster\'s mean binding energy deviates from a size-matched random distribution. Indicates statistical significance of coordinated binding patterns. Negative z-scores suggest functionally interacting protein complexes',
    'effect_size': 'Cohen\'s d-like metric quantifying the magnitude of binding energy effect. Measures practical significance of physical interactions independent of sample size. High absolute values indicate strong protein-protein interaction signals'
}

def train_model(df, features, feature_name=None, test_size=0.3, random_state=42):
    """Train classification model"""
    
    if feature_name:
        print(f"\n{'='*80}")
        print(f"INDIVIDUAL FEATURE MODEL: {feature_name}")
        print(f"Description: {FEATURE_DESCRIPTIONS[feature_name]}")
        print('='*80)
        features_to_use = [feature_name]
    else:
        print(f"\n{'='*80}")
        print(f"MULTI-FEATURE MODEL: {len(features)} features combined")
        print('='*80)
        features_to_use = features
    
    # Prepare data
    # Use 'name' if 'cluster_name' doesn't exist
    cluster_col = 'cluster_name' if 'cluster_name' in df.columns else 'name'
    data = df[[cluster_col] + features_to_use + ['label']].dropna()
    X = data[features_to_use].values
    y = data['label'].values
    
    print(f"\nData: {len(data)} samples ({y.sum()} MGC, {(y==0).sum()} Random)")
    
    # Train/Test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, stratify=y, random_state=random_state
    )
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train model
    model = LogisticRegression(random_state=random_state, max_iter=1000)
    model.fit(X_train_scaled, y_train)
    
    # Get weights
    weights = model.coef_[0]
    intercept = model.intercept_[0]
    
    print(f"\n📊 LEARNED WEIGHTS:")
    if len(features_to_use) == 1:
        print(f"  Weight: {weights[0]:.6f}")
        print(f"  Intercept: {intercept:.6f}")
    else:
        print(f"  Intercept: {intercept:.6f}")
        for feat, weight in zip(features_to_use, weights):
            print(f"  {feat:45s}: {weight:+.6f}")
        
        # Feature importance
        feature_importance = sorted(zip(features_to_use, np.abs(weights)), key=lambda x: x[1], reverse=True)
        print(f"\n🏆 FEATURE IMPORTANCE:")
        for rank, (feat, importance) in enumerate(feature_importance, 1):
            print(f"  {rank}. {feat:45s}: {importance:.6f}")
    
    # Predictions
    y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
    
    # ROC & PR
    fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
    roc_auc = auc(fpr, tpr)
    precision, recall, pr_thresholds = precision_recall_curve(y_test, y_pred_proba)
    pr_auc = auc(recall, precision)
    
    # Find F1-optimal threshold
    f1_scores = 2 * (precision[:-1] * recall[:-1]) / (precision[:-1] + recall[:-1] + 1e-10)
    f1_optimal_idx = np.argmax(f1_scores)
    f1_optimal_threshold = pr_thresholds[f1_optimal_idx]
    
    print(f"\n🎯 PERFORMANCE:")
    print(f"  ROC AUC: {roc_auc:.4f}")
    print(f"  PR AUC:  {pr_auc:.4f}")
    print(f"  F1-optimal threshold: {f1_optimal_threshold:.4f}")
    
    # Evaluate at F1-optimal threshold
    y_pred = (y_pred_proba >= f1_optimal_threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
    
    precision_val = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall_val = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = f1_score(y_test, y_pred)
    mcc = matthews_corrcoef(y_test, y_pred)
    
    print(f"\n📈 METRICS:")
    print(f"  Precision: {precision_val:.4f} | Recall: {recall_val:.4f} | F1: {f1:.4f} | MCC: {mcc:.4f}")
    print(f"  TP: {tp:4d} | FP: {fp:4d} | TN: {tn:5d} | FN: {fn:4d}")
    
    # Cross-validation
    cv_scores = cross_validate(model, scaler.transform(X), y, cv=5, 
                               scoring=['roc_auc', 'f1'], return_train_score=False)
    
    print(f"\n🔄 CROSS-VALIDATION:")
    print(f"  ROC AUC: {cv_scores['test_roc_auc'].mean():.4f} ± {cv_scores['test_roc_auc'].std():.4f}")
    print(f"  F1:      {cv_scores['test_f1'].mean():.4f} ± {cv_scores['test_f1'].std():.4f}")
    
    return {
        'features': features_to_use,
        'weights': weights,
        'intercept': intercept,
        'threshold': f1_optimal_threshold,
        'roc_auc': roc_auc,
        'pr_auc': pr_auc,
        'f1': f1,
        'mcc': mcc,
        'precision': precision_val,
        'recall': recall_val,
        'cv_roc_auc': cv_scores['test_roc_auc'].mean(),
        'cv_roc_auc_std': cv_scores['test_roc_auc'].std()
    }

def main():
    print("="*80)
    print("DOCKING ML CLASSIFICATION ANALYSIS")
    print("="*80)
    
    # Load data
    metrics_file = "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data/actual_random/lightdock_cluster_metrics.csv"
    df = pd.read_csv(metrics_file)
    # Create binary labels: 1 for BGC/MGC_CANDIDATE, 0 for RANDOM
    df['label'] = df['category'].apply(lambda x: 0 if x == 'RANDOM' else 1)
    
    print(f"\nLoaded {len(df)} clusters ({df['label'].sum()} MGC, {(df['label']==0).sum()} Random)")
    
    # Individual models
    print(f"\n{'='*80}")
    print("PART 1: INDIVIDUAL FEATURE MODELS")
    print("="*80)
    
    individual_results = []
    for feature in DOCKING_FEATURES:
        result = train_model(df, DOCKING_FEATURES, feature_name=feature)
        individual_results.append(result)
    
    # Multi-feature model
    print(f"\n{'='*80}")
    print("PART 2: MULTI-FEATURE MODEL")
    print("="*80)
    
    multi_result = train_model(df, DOCKING_FEATURES)
    
    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print("="*80)
    
    print(f"\nIndividual Features (ranked by ROC AUC):")
    sorted_results = sorted(individual_results, key=lambda x: x['roc_auc'], reverse=True)
    for i, r in enumerate(sorted_results, 1):
        print(f"  {i}. {r['features'][0]:45s} ROC: {r['roc_auc']:.4f} | F1: {r['f1']:.4f}")
    
    print(f"\nMulti-Feature Model:")
    print(f"  ROC AUC: {multi_result['roc_auc']:.4f} | F1: {multi_result['f1']:.4f} | MCC: {multi_result['mcc']:.4f}")
    
    # Save results
    summary = []
    for r in individual_results:
        summary.append({
            'model_type': 'individual',
            'feature': r['features'][0],
            'weight': r['weights'][0],
            'threshold': r['threshold'],
            'roc_auc': r['roc_auc'],
            'f1': r['f1'],
            'mcc': r['mcc']
        })
    
    summary.append({
        'model_type': 'multi_feature',
        'feature': 'ALL_COMBINED',
        'weight': None,
        'threshold': multi_result['threshold'],
        'roc_auc': multi_result['roc_auc'],
        'f1': multi_result['f1'],
        'mcc': multi_result['mcc']
    })
    
    summary_df = pd.DataFrame(summary)
    output_file = "/groups/itay_mayrose/alongonda/desktop/docking_ml_classification_summary.csv"
    summary_df.to_csv(output_file, index=False)
    print(f"\n💾 Saved: {output_file}")
    
    # Save multi-feature weights
    weights_df = pd.DataFrame({
        'feature': multi_result['features'],
        'weight': multi_result['weights'],
        'abs_importance': np.abs(multi_result['weights'])
    }).sort_values('abs_importance', ascending=False)
    
    weights_file = "/groups/itay_mayrose/alongonda/desktop/docking_multi_feature_weights.csv"
    weights_df.to_csv(weights_file, index=False)
    print(f"💾 Saved: {weights_file}")

if __name__ == "__main__":
    main()

