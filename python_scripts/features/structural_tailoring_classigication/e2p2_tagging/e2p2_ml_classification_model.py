#!/usr/bin/env python3
"""
E2P2 Feature Classification Models - Individual & Multi-feature
Uses proper ML frameworks (scikit-learn) for classification with learned weights
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import (roc_curve, auc, precision_recall_curve, confusion_matrix, 
                            classification_report, matthews_corrcoef, f1_score)
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score, cross_validate
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# E2P2 FEATURES
# ============================================================================
E2P2_FEATURES = [
    'num_distinct_enzyme_classes',
    'num_distinct_enzyme_subclasses',
    'num_distinct_enzyme_families',
    'num_distinct_enzyme_subfamilies',
    'total_ec_numbers'
]

FEATURE_DESCRIPTIONS = {
    'num_distinct_enzyme_classes': 'Count of distinct first-level EC numbers (format: X.-.-.-). Represents diversity across the 6 main enzyme classes: 1=Oxidoreductases, 2=Transferases, 3=Hydrolases, 4=Lyases, 5=Isomerases, 6=Ligases. Low values suggest pathway-specific enzymatic function',
    'num_distinct_enzyme_subclasses': 'Count of distinct second-level EC numbers (format: X.Y.-.-). Represents diversity of enzyme subclasses within main categories. Measures functional breadth at intermediate specificity',
    'num_distinct_enzyme_families': 'Count of distinct third-level EC numbers (format: X.Y.Z.-). Represents diversity of enzyme families with specific reaction types. Higher granularity than subclasses, indicates mechanistic diversity',
    'num_distinct_enzyme_subfamilies': 'Count of distinct fourth-level EC numbers (format: X.Y.Z.W). Most specific enzyme classification level, representing exact substrate-specific enzymes. High diversity suggests metabolically promiscuous clusters; low diversity suggests pathway specialization',
    'total_ec_numbers': 'Total count of all E2P2 enzyme predictions including duplicates. High counts with low diversity indicate gene duplications and pathway amplification typical of MGCs. Positive predictor: more enzymes = higher MGC likelihood'
}

def load_data(feature_file):
    """Load E2P2 feature data"""
    df = pd.read_csv(feature_file)
    print(f"Loaded {len(df)} clusters from file")
    
    # Create binary labels: 1 for MGC, 0 for RANDOM
    df['label'] = df['classification_label'].astype(int)
    
    print(f"  MGC clusters: {df['label'].sum()}")
    print(f"  Random clusters: {(df['label']==0).sum()}")
    
    return df

def train_individual_feature_model(df, feature_name, test_size=0.3, random_state=42):
    """
    Train classification model for a SINGLE feature
    Returns: trained model, predictions, metrics, and learned weight
    """
    print(f"\n{'='*80}")
    print(f"INDIVIDUAL FEATURE MODEL: {feature_name}")
    print(f"Description: {FEATURE_DESCRIPTIONS[feature_name]}")
    print('='*80)
    
    # Prepare data
    data = df[['cluster_name', feature_name, 'label']].dropna()
    X = data[[feature_name]].values
    y = data['label'].values
    
    print(f"\nData: {len(data)} samples ({y.sum()} MGC, {(y==0).sum()} Random)")
    
    # Descriptive statistics
    mgc_mean = X[y==1].mean()
    random_mean = X[y==0].mean()
    print(f"  MGC mean: {mgc_mean:.4f}")
    print(f"  Random mean: {random_mean:.4f}")
    print(f"  Direction: {'LOWER values = MGC' if mgc_mean < random_mean else 'HIGHER values = MGC'}")
    
    # Train/Test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, stratify=y, random_state=random_state
    )
    
    print(f"\nTrain/Test Split:")
    print(f"  Train: {len(X_train)} samples ({y_train.sum()} MGC, {(y_train==0).sum()} Random)")
    print(f"  Test:  {len(X_test)} samples ({y_test.sum()} MGC, {(y_test==0).sum()} Random)")
    
    # Train Logistic Regression model
    model = LogisticRegression(random_state=random_state, max_iter=1000)
    model.fit(X_train, y_train)
    
    # Get learned weight (coefficient)
    weight = model.coef_[0][0]
    intercept = model.intercept_[0]
    
    print(f"\n📊 LEARNED MODEL PARAMETERS:")
    print(f"  Weight (coefficient): {weight:.6f}")
    print(f"  Intercept: {intercept:.6f}")
    print(f"  Formula: P(MGC) = sigmoid({weight:.6f} * {feature_name} + {intercept:.6f})")
    
    # Predict probabilities on test set
    y_pred_proba = model.predict_proba(X_test)[:, 1]  # Probability of being MGC
    
    # ROC Analysis
    fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
    roc_auc = auc(fpr, tpr)
    
    # Precision-Recall Analysis
    precision, recall, pr_thresholds = precision_recall_curve(y_test, y_pred_proba)
    pr_auc = auc(recall, precision)
    
    # Find optimal threshold using Youden's J
    youden_j = tpr - fpr
    optimal_idx = np.argmax(youden_j)
    optimal_threshold = thresholds[optimal_idx]
    
    # Find F1-optimal threshold
    f1_scores = 2 * (precision[:-1] * recall[:-1]) / (precision[:-1] + recall[:-1] + 1e-10)
    f1_optimal_idx = np.argmax(f1_scores)
    f1_optimal_threshold = pr_thresholds[f1_optimal_idx]
    
    print(f"\n🎯 TEST SET PERFORMANCE:")
    print(f"  ROC AUC: {roc_auc:.4f}")
    print(f"  PR AUC:  {pr_auc:.4f}")
    
    print(f"\n🔧 OPTIMAL CLASSIFICATION THRESHOLDS:")
    print(f"  Youden's J optimal: {optimal_threshold:.4f}")
    print(f"  F1-optimal: {f1_optimal_threshold:.4f}")
    
    # Evaluate at F1-optimal threshold
    y_pred = (y_pred_proba >= f1_optimal_threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
    
    precision_val = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall_val = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    f1 = f1_score(y_test, y_pred)
    mcc = matthews_corrcoef(y_test, y_pred)
    accuracy = (tp + tn) / (tp + tn + fp + fn)
    
    print(f"\n📈 CLASSIFICATION METRICS (F1-optimal threshold):")
    print(f"  Accuracy:    {accuracy:.4f}")
    print(f"  Precision:   {precision_val:.4f}")
    print(f"  Recall:      {recall_val:.4f}")
    print(f"  Specificity: {specificity:.4f}")
    print(f"  F1-Score:    {f1:.4f}")
    print(f"  MCC:         {mcc:.4f}")
    
    print(f"\n📋 CONFUSION MATRIX:")
    print(f"  True Negatives:  {tn:5d}  |  False Positives: {fp:5d}")
    print(f"  False Negatives: {fn:5d}  |  True Positives:  {tp:5d}")
    
    # Cross-validation
    print(f"\n🔄 5-FOLD CROSS-VALIDATION:")
    cv_scores = cross_validate(
        model, X, y, cv=5,
        scoring=['roc_auc', 'precision', 'recall', 'f1'],
        return_train_score=False
    )
    
    print(f"  ROC AUC:   {cv_scores['test_roc_auc'].mean():.4f} ± {cv_scores['test_roc_auc'].std():.4f}")
    print(f"  Precision: {cv_scores['test_precision'].mean():.4f} ± {cv_scores['test_precision'].std():.4f}")
    print(f"  Recall:    {cv_scores['test_recall'].mean():.4f} ± {cv_scores['test_recall'].std():.4f}")
    print(f"  F1-Score:  {cv_scores['test_f1'].mean():.4f} ± {cv_scores['test_f1'].std():.4f}")
    
    # Production recommendation
    if roc_auc >= 0.9 and cv_scores['test_roc_auc'].std() < 0.05:
        recommendation = "✅ EXCELLENT - Highly recommended for production"
    elif roc_auc >= 0.8 and cv_scores['test_roc_auc'].std() < 0.1:
        recommendation = "✓ GOOD - Suitable for production with monitoring"
    elif roc_auc >= 0.7:
        recommendation = "⚠ FAIR - Use with caution"
    else:
        recommendation = "❌ POOR - Not recommended for production"
    
    print(f"\n💡 RECOMMENDATION: {recommendation}")
    
    return {
        'feature_name': feature_name,
        'model': model,
        'weight': weight,
        'intercept': intercept,
        'scaler': None,
        'optimal_threshold_youden': optimal_threshold,
        'optimal_threshold_f1': f1_optimal_threshold,
        'test_roc_auc': roc_auc,
        'test_pr_auc': pr_auc,
        'test_accuracy': accuracy,
        'test_precision': precision_val,
        'test_recall': recall_val,
        'test_f1': f1,
        'test_mcc': mcc,
        'cv_roc_auc_mean': cv_scores['test_roc_auc'].mean(),
        'cv_roc_auc_std': cv_scores['test_roc_auc'].std(),
        'cv_f1_mean': cv_scores['test_f1'].mean(),
        'cv_f1_std': cv_scores['test_f1'].std(),
        'X_test': X_test,
        'y_test': y_test,
        'y_pred_proba': y_pred_proba
    }

def train_multi_feature_model(df, features=None, test_size=0.3, random_state=42):
    """
    Train classification model using MULTIPLE features
    Returns: trained model, predictions, metrics, and learned weights for each feature
    """
    if features is None:
        features = E2P2_FEATURES
    
    print(f"\n{'='*80}")
    print(f"MULTI-FEATURE CLASSIFICATION MODEL")
    print(f"Using {len(features)} features combined")
    print('='*80)
    
    for i, feat in enumerate(features, 1):
        print(f"  {i}. {feat}: {FEATURE_DESCRIPTIONS[feat]}")
    
    # Prepare data
    data = df[['cluster_name'] + features + ['label']].dropna()
    X = data[features].values
    y = data['label'].values
    
    print(f"\nData: {len(data)} samples ({y.sum()} MGC, {(y==0).sum()} Random)")
    
    # Train/Test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, stratify=y, random_state=random_state
    )
    
    print(f"\nTrain/Test Split:")
    print(f"  Train: {len(X_train)} samples ({y_train.sum()} MGC, {(y_train==0).sum()} Random)")
    print(f"  Test:  {len(X_test)} samples ({y_test.sum()} MGC, {(y_test==0).sum()} Random)")
    
    # Feature scaling for better interpretability
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train Logistic Regression model
    model = LogisticRegression(random_state=random_state, max_iter=1000)
    model.fit(X_train_scaled, y_train)
    
    # Get learned weights (coefficients)
    weights = model.coef_[0]
    intercept = model.intercept_[0]
    
    print(f"\n📊 LEARNED MODEL PARAMETERS (on standardized features):")
    print(f"  Intercept: {intercept:.6f}")
    print(f"\n  Feature Weights (coefficients):")
    for feat, weight in zip(features, weights):
        direction = "↑ Positive" if weight > 0 else "↓ Negative"
        importance = abs(weight)
        print(f"    {feat:45s}: {weight:+.6f} ({direction}, importance: {importance:.6f})")
    
    # Feature importance ranking
    feature_importance = sorted(zip(features, np.abs(weights)), key=lambda x: x[1], reverse=True)
    print(f"\n  Feature Importance Ranking:")
    for rank, (feat, importance) in enumerate(feature_importance, 1):
        print(f"    {rank}. {feat:45s}: {importance:.6f}")
    
    # Predict probabilities on test set
    y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
    
    # ROC Analysis
    fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
    roc_auc = auc(fpr, tpr)
    
    # Precision-Recall Analysis
    precision, recall, pr_thresholds = precision_recall_curve(y_test, y_pred_proba)
    pr_auc = auc(recall, precision)
    
    # Find optimal threshold using Youden's J
    youden_j = tpr - fpr
    optimal_idx = np.argmax(youden_j)
    optimal_threshold = thresholds[optimal_idx]
    
    # Find F1-optimal threshold
    f1_scores = 2 * (precision[:-1] * recall[:-1]) / (precision[:-1] + recall[:-1] + 1e-10)
    f1_optimal_idx = np.argmax(f1_scores)
    f1_optimal_threshold = pr_thresholds[f1_optimal_idx]
    
    print(f"\n🎯 TEST SET PERFORMANCE:")
    print(f"  ROC AUC: {roc_auc:.4f}")
    print(f"  PR AUC:  {pr_auc:.4f}")
    
    print(f"\n🔧 OPTIMAL CLASSIFICATION THRESHOLDS:")
    print(f"  Youden's J optimal: {optimal_threshold:.4f}")
    print(f"  F1-optimal: {f1_optimal_threshold:.4f}")
    
    # Evaluate at F1-optimal threshold
    y_pred = (y_pred_proba >= f1_optimal_threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
    
    precision_val = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall_val = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    f1 = f1_score(y_test, y_pred)
    mcc = matthews_corrcoef(y_test, y_pred)
    accuracy = (tp + tn) / (tp + tn + fp + fn)
    
    print(f"\n📈 CLASSIFICATION METRICS (F1-optimal threshold):")
    print(f"  Accuracy:    {accuracy:.4f}")
    print(f"  Precision:   {precision_val:.4f}")
    print(f"  Recall:      {recall_val:.4f}")
    print(f"  Specificity: {specificity:.4f}")
    print(f"  F1-Score:    {f1:.4f}")
    print(f"  MCC:         {mcc:.4f}")
    
    print(f"\n📋 CONFUSION MATRIX:")
    print(f"  True Negatives:  {tn:5d}  |  False Positives: {fp:5d}")
    print(f"  False Negatives: {fn:5d}  |  True Positives:  {tp:5d}")
    
    # Cross-validation
    print(f"\n🔄 5-FOLD CROSS-VALIDATION:")
    cv_scores = cross_validate(
        model, scaler.transform(X), y, cv=5,
        scoring=['roc_auc', 'precision', 'recall', 'f1'],
        return_train_score=False
    )
    
    print(f"  ROC AUC:   {cv_scores['test_roc_auc'].mean():.4f} ± {cv_scores['test_roc_auc'].std():.4f}")
    print(f"  Precision: {cv_scores['test_precision'].mean():.4f} ± {cv_scores['test_precision'].std():.4f}")
    print(f"  Recall:    {cv_scores['test_recall'].mean():.4f} ± {cv_scores['test_recall'].std():.4f}")
    print(f"  F1-Score:  {cv_scores['test_f1'].mean():.4f} ± {cv_scores['test_f1'].std():.4f}")
    
    # Production recommendation
    if roc_auc >= 0.9 and cv_scores['test_roc_auc'].std() < 0.05:
        recommendation = "✅ EXCELLENT - Highly recommended for production"
    elif roc_auc >= 0.8 and cv_scores['test_roc_auc'].std() < 0.1:
        recommendation = "✓ GOOD - Suitable for production with monitoring"
    elif roc_auc >= 0.7:
        recommendation = "⚠ FAIR - Use with caution"
    else:
        recommendation = "❌ POOR - Not recommended for production"
    
    print(f"\n💡 RECOMMENDATION: {recommendation}")
    
    return {
        'model_type': 'multi_feature',
        'features': features,
        'model': model,
        'scaler': scaler,
        'weights': weights,
        'intercept': intercept,
        'feature_importance': feature_importance,
        'optimal_threshold_youden': optimal_threshold,
        'optimal_threshold_f1': f1_optimal_threshold,
        'test_roc_auc': roc_auc,
        'test_pr_auc': pr_auc,
        'test_accuracy': accuracy,
        'test_precision': precision_val,
        'test_recall': recall_val,
        'test_f1': f1,
        'test_mcc': mcc,
        'cv_roc_auc_mean': cv_scores['test_roc_auc'].mean(),
        'cv_roc_auc_std': cv_scores['test_roc_auc'].std(),
        'cv_f1_mean': cv_scores['test_f1'].mean(),
        'cv_f1_std': cv_scores['test_f1'].std(),
        'X_test': X_test_scaled,
        'y_test': y_test,
        'y_pred_proba': y_pred_proba
    }

def main():
    """Main analysis function"""
    
    print("="*80)
    print("E2P2 CLASSIFICATION MODEL ANALYSIS")
    print("Individual & Multi-Feature ML Models with Learned Weights")
    print("="*80)
    
    # Load data
    feature_file = "/groups/itay_mayrose/alongonda/Plant_MGC/fixed_kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/random_kegg_e2p2_feature_extraction_results/e2p2_feature_extraction_detailed.csv"
    
    df = load_data(feature_file)
    
    # ========================================================================
    # PART 1: INDIVIDUAL FEATURE MODELS
    # ========================================================================
    print("\n" + "="*80)
    print("PART 1: INDIVIDUAL FEATURE CLASSIFICATION MODELS")
    print("="*80)
    
    individual_results = {}
    for feature in E2P2_FEATURES:
        result = train_individual_feature_model(df, feature)
        individual_results[feature] = result
    
    # ========================================================================
    # PART 2: MULTI-FEATURE MODEL
    # ========================================================================
    print("\n" + "="*80)
    print("PART 2: MULTI-FEATURE CLASSIFICATION MODEL")
    print("="*80)
    
    multi_result = train_multi_feature_model(df, E2P2_FEATURES)
    
    # ========================================================================
    # FINAL COMPARISON
    # ========================================================================
    print("\n" + "="*80)
    print("FINAL COMPARISON: ALL MODELS")
    print("="*80)
    
    print("\n📊 Individual Feature Models - Ranked by Test ROC AUC:")
    print("-" * 80)
    
    # Sort by ROC AUC
    sorted_results = sorted(individual_results.items(), 
                          key=lambda x: x[1]['test_roc_auc'], 
                          reverse=True)
    
    for rank, (feature, result) in enumerate(sorted_results, 1):
        print(f"\n{rank}. {feature}")
        print(f"   Weight: {result['weight']:+.6f} | Threshold: {result['optimal_threshold_f1']:.4f}")
        print(f"   ROC AUC: {result['test_roc_auc']:.4f} | F1: {result['test_f1']:.4f} | MCC: {result['test_mcc']:.4f}")
        print(f"   CV Stability: {result['cv_roc_auc_std']:.4f}")
    
    print(f"\n{'='*80}")
    print("Multi-Feature Combined Model:")
    print(f"   ROC AUC: {multi_result['test_roc_auc']:.4f} | F1: {multi_result['test_f1']:.4f} | MCC: {multi_result['test_mcc']:.4f}")
    print(f"   Threshold: {multi_result['optimal_threshold_f1']:.4f}")
    print(f"   CV Stability: {multi_result['cv_roc_auc_std']:.4f}")
    
    # Save summary results
    summary_data = []
    
    for feature, result in individual_results.items():
        summary_data.append({
            'model_type': 'individual',
            'feature': feature,
            'weight': result['weight'],
            'intercept': result['intercept'],
            'optimal_threshold': result['optimal_threshold_f1'],
            'test_roc_auc': result['test_roc_auc'],
            'test_pr_auc': result['test_pr_auc'],
            'test_f1': result['test_f1'],
            'test_mcc': result['test_mcc'],
            'cv_roc_auc_mean': result['cv_roc_auc_mean'],
            'cv_roc_auc_std': result['cv_roc_auc_std']
        })
    
    # Add multi-feature model
    summary_data.append({
        'model_type': 'multi_feature',
        'feature': 'ALL_FEATURES_COMBINED',
        'weight': None,
        'intercept': multi_result['intercept'],
        'optimal_threshold': multi_result['optimal_threshold_f1'],
        'test_roc_auc': multi_result['test_roc_auc'],
        'test_pr_auc': multi_result['test_pr_auc'],
        'test_f1': multi_result['test_f1'],
        'test_mcc': multi_result['test_mcc'],
        'cv_roc_auc_mean': multi_result['cv_roc_auc_mean'],
        'cv_roc_auc_std': multi_result['cv_roc_auc_std']
    })
    
    summary_df = pd.DataFrame(summary_data)
    output_file = "/groups/itay_mayrose/alongonda/desktop/random_kegg_e2p2_classification_models_summary.csv"
    summary_df.to_csv(output_file, index=False)
    
    print(f"\n💾 Summary results saved to: {output_file}")
    
    # Save multi-feature weights separately
    weights_df = pd.DataFrame({
        'feature': multi_result['features'],
        'weight': multi_result['weights'],
        'abs_importance': np.abs(multi_result['weights'])
    }).sort_values('abs_importance', ascending=False)
    
    weights_output = "/groups/itay_mayrose/alongonda/desktop/random_kegg_e2p2_multi_feature_weights.csv"
    weights_df.to_csv(weights_output, index=False)
    print(f"💾 Multi-feature weights saved to: {weights_output}")
    
    return individual_results, multi_result

if __name__ == "__main__":
    individual_results, multi_result = main()

