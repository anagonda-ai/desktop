#!/usr/bin/env python3
"""
Comprehensive ML Classification Analysis - Incremental Feature Addition
Uses the same ML approach as individual models: Logistic Regression with learned weights
"""

import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (roc_curve, auc, precision_recall_curve, confusion_matrix, 
                            matthews_corrcoef, f1_score)
from sklearn.model_selection import train_test_split, cross_validate
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

class ComprehensiveMLAnalyzer:
    def __init__(self):
        BASE_DIR = "/groups/itay_mayrose/alongonda/desktop/python_scripts/features/final_data"
        
        # Feature configurations (cleaned - no approximations)
        self.feature_configs = {
            'docking': {
                'path': f'{BASE_DIR}/lightdock_cluster_metrics.csv',
                'cluster_col': 'name',
                'features': [
                    'mean_score_non_self',
                    'enrichment_score',
                    'z_score',
                    'effect_size'
                ],
                'name': 'Docking'
            },
            'foldseek': {
                'path': f'{BASE_DIR}/foldseek_cluster_metrics_enhanced.csv',
                'cluster_col': 'name',
                'features': [
                    'mean_score_non_self',
                    'enrichment_score',
                    'z_score',
                    'effect_size',
                    'foldseek_match_coverage'
                ],
                'name': 'Foldseek'
            },
            'e2p2': {
                'path': f'{BASE_DIR}/e2p2_feature_extraction_detailed.csv',
                'cluster_col': 'cluster_name',
                'features': [
                    'num_distinct_enzyme_classes',
                    'num_distinct_enzyme_subfamilies',
                    'total_ec_numbers'
                ],
                'name': 'E2P2'
            },
            'cladepp': {
                'path': f'{BASE_DIR}/cladepp_cluster_metrics_enhanced.csv',
                'cluster_col': 'cluster_name',
                'features': [
                    'mean_cladepp_score',
                    'weighted_cladepp_score',
                    'positive_correlation_ratio',
                    'cladepp_multi_clade_high',
                    'cladepp_multi_clade_medium',
                    'cladepp_conservation_consistency',
                    'cladepp_max_pair_score'
                ],
                'name': 'Cladepp'
            }
        }
        
        print("="*90)
        print("COMPREHENSIVE ML CLASSIFICATION ANALYSIS")
        print("Incremental Feature Addition with Learned Weights")
        print("="*90)
        
        total_features = sum(len(c['features']) for c in self.feature_configs.values())
        print(f"\nTotal features across all classes: {total_features}")
        for name, config in self.feature_configs.items():
            print(f"  {config['name']:12s}: {len(config['features']):2d} features")
    
    def load_and_merge_data(self):
        """Load all feature data and merge on cluster names"""
        print("\n" + "="*90)
        print("LOADING DATA")
        print("="*90)
        
        merged_df = None
        
        for key, config in self.feature_configs.items():
            df = pd.read_csv(config['path'])
            
            # Create label
            if 'category' in df.columns:
                df['label'] = df['category'].apply(lambda x: 0 if x == 'RANDOM' else 1)
            elif 'classification_label' in df.columns:
                df['label'] = df['classification_label'].astype(int)
            else:
                raise ValueError(f"No label column found in {key}")
            
            # Select columns
            cluster_col = config['cluster_col']
            cols_to_keep = [cluster_col, 'label'] + config['features']
            df = df[cols_to_keep]
            
            # Rename cluster column to standard name
            df = df.rename(columns={cluster_col: 'cluster_name'})
            
            # Add prefix to feature names to avoid collisions
            prefix = key + '_'
            rename_dict = {feat: prefix + feat for feat in config['features']}
            df = df.rename(columns=rename_dict)
            
            # Update config with prefixed names
            config['prefixed_features'] = [prefix + feat for feat in config['features']]
            
            print(f"\n{config['name']:12s}: {len(df):6d} clusters, {len(config['features']):2d} features")
            
            # Merge
            if merged_df is None:
                merged_df = df
            else:
                merged_df = merged_df.merge(df, on=['cluster_name', 'label'], how='inner')
        
        print(f"\nMerged dataset: {len(merged_df)} clusters")
        print(f"  MGC:    {merged_df['label'].sum()}")
        print(f"  Random: {(merged_df['label']==0).sum()}")
        
        return merged_df
    
    def train_model(self, X, y, feature_names, test_size=0.3, random_state=42):
        """Train Logistic Regression model"""
        # Train/Test split
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, stratify=y, random_state=random_state
        )
        
        # Scale
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # Train
        model = LogisticRegression(random_state=random_state, max_iter=1000)
        model.fit(X_train_scaled, y_train)
        
        # Predict
        y_pred_proba = model.predict_proba(X_test_scaled)[:, 1]
        
        # ROC & PR
        fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
        roc_auc = auc(fpr, tpr)
        precision, recall, pr_thresholds = precision_recall_curve(y_test, y_pred_proba)
        pr_auc = auc(recall, precision)
        
        # F1-optimal threshold
        f1_scores = 2 * (precision[:-1] * recall[:-1]) / (precision[:-1] + recall[:-1] + 1e-10)
        f1_optimal_idx = np.argmax(f1_scores)
        f1_optimal_threshold = pr_thresholds[f1_optimal_idx]
        
        # Evaluate
        y_pred = (y_pred_proba >= f1_optimal_threshold).astype(int)
        tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
        
        precision_val = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall_val = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = f1_score(y_test, y_pred)
        mcc = matthews_corrcoef(y_test, y_pred)
        
        # Cross-validation
        cv_scores = cross_validate(model, scaler.transform(X), y, cv=5, 
                                   scoring=['roc_auc', 'f1'], return_train_score=False)
        
        return {
            'model': model,
            'scaler': scaler,
            'weights': model.coef_[0],
            'intercept': model.intercept_[0],
            'feature_names': feature_names,
            'threshold': f1_optimal_threshold,
            'roc_auc': roc_auc,
            'pr_auc': pr_auc,
            'f1': f1,
            'mcc': mcc,
            'precision': precision_val,
            'recall': recall_val,
            'cv_roc_auc_mean': cv_scores['test_roc_auc'].mean(),
            'cv_roc_auc_std': cv_scores['test_roc_auc'].std(),
            'cv_f1_mean': cv_scores['test_f1'].mean(),
            'cv_f1_std': cv_scores['test_f1'].std()
        }
    
    def incremental_analysis(self, df):
        """Perform incremental feature addition analysis"""
        print("\n" + "="*90)
        print("INCREMENTAL FEATURE ADDITION ANALYSIS")
        print("="*90)
        
        # Get all features in order of feature class
        all_features = []
        feature_class_map = {}
        
        for key, config in self.feature_configs.items():
            for feat in config['prefixed_features']:
                all_features.append(feat)
                feature_class_map[feat] = config['name']
        
        y = df['label'].values
        results = []
        
        # Train models with incremental features
        for i in range(1, len(all_features) + 1):
            features_subset = all_features[:i]
            X = df[features_subset].values
            
            # Remove rows with NaN
            valid_mask = ~np.isnan(X).any(axis=1)
            X_clean = X[valid_mask]
            y_clean = y[valid_mask]
            
            if len(np.unique(y_clean)) < 2:
                continue
            
            result = self.train_model(X_clean, y_clean, features_subset)
            
            results.append({
                'n_features': i,
                'features': ', '.join(features_subset) if i <= 3 else f"{features_subset[0]}, ... [{i} total]",
                'last_added': features_subset[-1],
                'feature_class': feature_class_map[features_subset[-1]],
                'roc_auc': result['roc_auc'],
                'f1': result['f1'],
                'mcc': result['mcc'],
                'precision': result['precision'],
                'recall': result['recall'],
                'cv_roc_auc': result['cv_roc_auc_mean'],
                'cv_roc_std': result['cv_roc_auc_std']
            })
            
            print(f"\n{i:2d} features | Last added: {features_subset[-1]:45s} ({feature_class_map[features_subset[-1]]:10s})")
            print(f"    ROC AUC: {result['roc_auc']:.4f} | F1: {result['f1']:.4f} | MCC: {result['mcc']:.4f}")
        
        return pd.DataFrame(results)
    
    def train_individual_class_models(self, df):
        """Train individual models for each feature class"""
        print("\n" + "="*90)
        print("INDIVIDUAL FEATURE CLASS MODELS")
        print("="*90)
        
        class_results = []
        
        for key, config in self.feature_configs.items():
            print(f"\n{'='*80}")
            print(f"{config['name'].upper()} MODEL ({len(config['features'])} features)")
            print('='*80)
            
            features = config['prefixed_features']
            data = df[['cluster_name'] + features + ['label']].dropna()
            X = data[features].values
            y = data['label'].values
            
            print(f"Data: {len(data)} samples ({y.sum()} MGC, {(y==0).sum()} Random)")
            
            result = self.train_model(X, y, features)
            
            print(f"\n📊 WEIGHTS:")
            print(f"  Intercept: {result['intercept']:.4f}")
            for feat, weight in zip(features, result['weights']):
                print(f"  {feat:45s}: {weight:+.4f}")
            
            # Feature importance
            importance = sorted(zip(features, np.abs(result['weights'])), key=lambda x: x[1], reverse=True)
            print(f"\n🏆 FEATURE IMPORTANCE:")
            for rank, (feat, imp) in enumerate(importance, 1):
                print(f"  {rank}. {feat:45s}: {imp:.4f}")
            
            print(f"\n🎯 PERFORMANCE:")
            print(f"  ROC AUC: {result['roc_auc']:.4f} | F1: {result['f1']:.4f} | MCC: {result['mcc']:.4f}")
            print(f"  Precision: {result['precision']:.4f} | Recall: {result['recall']:.4f}")
            print(f"  CV ROC AUC: {result['cv_roc_auc_mean']:.4f} ± {result['cv_roc_auc_std']:.4f}")
            
            class_results.append({
                'feature_class': config['name'],
                'n_features': len(features),
                'roc_auc': result['roc_auc'],
                'f1': result['f1'],
                'mcc': result['mcc'],
                'cv_roc_auc': result['cv_roc_auc_mean'],
                'top_feature': importance[0][0],
                'top_weight': importance[0][1]
            })
        
        return pd.DataFrame(class_results)
    
    def train_combined_model(self, df):
        """Train model with ALL features combined"""
        print("\n" + "="*90)
        print("COMBINED MODEL - ALL FEATURES")
        print("="*90)
        
        all_features = []
        for config in self.feature_configs.values():
            all_features.extend(config['prefixed_features'])
        
        data = df[['cluster_name'] + all_features + ['label']].dropna()
        X = data[all_features].values
        y = data['label'].values
        
        print(f"\nData: {len(data)} samples ({y.sum()} MGC, {(y==0).sum()} Random)")
        print(f"Features: {len(all_features)}")
        
        result = self.train_model(X, y, all_features)
        
        print(f"\n📊 WEIGHTS (Top 10):")
        print(f"  Intercept: {result['intercept']:.4f}")
        importance = sorted(zip(all_features, np.abs(result['weights'])), key=lambda x: x[1], reverse=True)
        for rank, (feat, imp) in enumerate(importance[:10], 1):
            weight = result['weights'][all_features.index(feat)]
            print(f"  {rank:2d}. {feat:45s}: {weight:+.4f} (importance: {imp:.4f})")
        
        print(f"\n🎯 PERFORMANCE:")
        print(f"  ROC AUC: {result['roc_auc']:.4f}")
        print(f"  PR AUC:  {result['pr_auc']:.4f}")
        print(f"  F1:      {result['f1']:.4f}")
        print(f"  MCC:     {result['mcc']:.4f}")
        print(f"  Precision: {result['precision']:.4f}")
        print(f"  Recall:    {result['recall']:.4f}")
        print(f"\n  CV ROC AUC: {result['cv_roc_auc_mean']:.4f} ± {result['cv_roc_auc_std']:.4f}")
        print(f"  CV F1:      {result['cv_f1_mean']:.4f} ± {result['cv_f1_std']:.4f}")
        
        return result, importance
    
    def run_full_analysis(self):
        """Run complete analysis pipeline"""
        # Load data
        df = self.load_and_merge_data()
        
        # Individual class models
        class_results = self.train_individual_class_models(df)
        
        # Incremental analysis
        incremental_results = self.incremental_analysis(df)
        
        # Combined model
        combined_result, feature_importance = self.train_combined_model(df)
        
        # Summary
        print("\n" + "="*90)
        print("FINAL SUMMARY")
        print("="*90)
        
        print("\n📊 Individual Feature Class Performance:")
        print("-" * 90)
        class_results_sorted = class_results.sort_values('roc_auc', ascending=False)
        for _, row in class_results_sorted.iterrows():
            print(f"  {row['feature_class']:12s} ({row['n_features']:2d} features): "
                  f"ROC AUC: {row['roc_auc']:.4f} | F1: {row['f1']:.4f} | "
                  f"Top: {row['top_feature'][:35]:35s}")
        
        print(f"\n📊 Combined Model Performance:")
        print(f"  ALL FEATURES: ROC AUC: {combined_result['roc_auc']:.4f} | "
              f"F1: {combined_result['f1']:.4f} | MCC: {combined_result['mcc']:.4f}")
        
        # Save results
        class_results.to_csv("/groups/itay_mayrose/alongonda/desktop/ml_individual_class_results.csv", index=False)
        incremental_results.to_csv("/groups/itay_mayrose/alongonda/desktop/ml_incremental_results.csv", index=False)
        
        importance_df = pd.DataFrame(feature_importance, columns=['feature', 'importance'])
        importance_df.to_csv("/groups/itay_mayrose/alongonda/desktop/ml_combined_feature_importance.csv", index=False)
        
        print("\n💾 Results saved:")
        print("  - ml_individual_class_results.csv")
        print("  - ml_incremental_results.csv")
        print("  - ml_combined_feature_importance.csv")
        
        return {
            'class_results': class_results,
            'incremental_results': incremental_results,
            'combined_result': combined_result,
            'feature_importance': feature_importance
        }

def main():
    analyzer = ComprehensiveMLAnalyzer()
    results = analyzer.run_full_analysis()
    return results

if __name__ == "__main__":
    results = main()
