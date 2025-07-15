"""
Advanced statistical analysis framework for bioinformatics.

This module provides comprehensive statistical tools including multiple testing
correction, enrichment analysis, and specialized bioinformatics statistics.
"""

import numpy as np
import pandas as pd
from typing import Any, Dict, List, Optional, Tuple, Union
from dataclasses import dataclass
from enum import Enum
import scipy.stats as stats
from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.contingency_tables import mcnemar
import networkx as nx
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import warnings
import logging

from ..core.exceptions import ValidationError, ProcessingError
from ..core.types import AnalysisResult, ProcessingStatus


class CorrectionMethod(Enum):
    """Multiple testing correction methods."""
    BONFERRONI = "bonferroni"
    SIDAK = "sidak" 
    HOLM_SIDAK = "holm-sidak"
    HOLM = "holm"
    SIMES_HOCHBERG = "simes-hochberg"
    HOMMEL = "hommel"
    FDR_BH = "fdr_bh"  # Benjamini-Hochberg
    FDR_BY = "fdr_by"  # Benjamini-Yekutieli
    FDR_TWO_STAGE = "fdr_tsbh"


class StatisticalTest(Enum):
    """Supported statistical tests."""
    T_TEST = "t_test"
    WILCOXON = "wilcoxon"
    MANN_WHITNEY = "mann_whitney"
    CHI_SQUARE = "chi_square"
    FISHER_EXACT = "fisher_exact"
    BINOMIAL = "binomial"
    HYPERGEOMETRIC = "hypergeometric"
    KOLMOGOROV_SMIRNOV = "ks_test"
    ANDERSON_DARLING = "anderson_darling"


@dataclass
class StatisticalResult:
    """Result of a statistical test."""
    test_name: str
    statistic: float
    p_value: float
    effect_size: Optional[float] = None
    confidence_interval: Optional[Tuple[float, float]] = None
    sample_size: Optional[int] = None
    degrees_freedom: Optional[int] = None
    interpretation: Optional[str] = None
    
    def is_significant(self, alpha: float = 0.05) -> bool:
        """Check if result is statistically significant."""
        return self.p_value < alpha


@dataclass
class EnrichmentResult:
    """Result of enrichment analysis."""
    term: str
    observed: int
    expected: float
    fold_enrichment: float
    p_value: float
    adjusted_p_value: Optional[float] = None
    genes: Optional[List[str]] = None
    
    def is_significant(self, alpha: float = 0.05, use_adjusted: bool = True) -> bool:
        """Check if enrichment is significant."""
        p_val = self.adjusted_p_value if use_adjusted and self.adjusted_p_value else self.p_value
        return p_val < alpha


class MultipleTestingCorrector:
    """Handle multiple testing corrections with various methods."""
    
    def __init__(self, method: CorrectionMethod = CorrectionMethod.FDR_BH):
        self.method = method
        self.logger = logging.getLogger(__name__)
    
    def correct_pvalues(
        self, 
        p_values: List[float], 
        alpha: float = 0.05
    ) -> Tuple[List[bool], List[float]]:
        """
        Apply multiple testing correction to p-values.
        
        Args:
            p_values: List of p-values to correct
            alpha: Significance threshold
            
        Returns:
            Tuple of (rejected hypotheses, corrected p-values)
        """
        if not p_values:
            return [], []
        
        try:
            rejected, corrected_pvals, _, _ = multipletests(
                p_values,
                alpha=alpha,
                method=self.method.value
            )
            
            self.logger.debug(
                f"Applied {self.method.value} correction to {len(p_values)} p-values. "
                f"Rejected {sum(rejected)} hypotheses at α={alpha}"
            )
            
            return rejected.tolist(), corrected_pvals.tolist()
            
        except Exception as e:
            raise ProcessingError(
                f"Multiple testing correction failed: {e}",
                processor_name=self.__class__.__name__,
                stage="correction"
            )
    
    def correct_dataframe(
        self, 
        df: pd.DataFrame, 
        pvalue_col: str, 
        alpha: float = 0.05,
        rejected_col: str = "significant",
        corrected_col: str = "adjusted_pvalue"
    ) -> pd.DataFrame:
        """
        Apply correction to DataFrame with p-values.
        
        Args:
            df: DataFrame containing p-values
            pvalue_col: Name of p-value column
            alpha: Significance threshold
            rejected_col: Name for significant results column
            corrected_col: Name for corrected p-values column
            
        Returns:
            DataFrame with correction results added
        """
        result_df = df.copy()
        
        if pvalue_col not in df.columns:
            raise ValidationError(f"P-value column '{pvalue_col}' not found in DataFrame")
        
        p_values = df[pvalue_col].fillna(1.0).tolist()  # Fill NaN with 1.0
        rejected, corrected = self.correct_pvalues(p_values, alpha)
        
        result_df[rejected_col] = rejected
        result_df[corrected_col] = corrected
        
        return result_df


class EnrichmentAnalyzer:
    """Perform enrichment analysis for gene sets and pathways."""
    
    def __init__(self, correction_method: CorrectionMethod = CorrectionMethod.FDR_BH):
        self.corrector = MultipleTestingCorrector(correction_method)
        self.logger = logging.getLogger(__name__)
    
    def hypergeometric_enrichment(
        self,
        gene_set: List[str],
        pathway_genes: Dict[str, List[str]],
        total_genes: int,
        min_overlap: int = 3
    ) -> List[EnrichmentResult]:
        """
        Perform hypergeometric enrichment analysis.
        
        Args:
            gene_set: List of genes of interest
            pathway_genes: Dictionary mapping pathway names to gene lists
            total_genes: Total number of genes in background
            min_overlap: Minimum overlap for testing
            
        Returns:
            List of enrichment results
        """
        results = []
        gene_set = set(gene_set)
        
        for pathway, pathway_gene_list in pathway_genes.items():
            pathway_set = set(pathway_gene_list)
            overlap = gene_set.intersection(pathway_set)
            
            if len(overlap) < min_overlap:
                continue
            
            # Hypergeometric test parameters
            k = len(overlap)  # Observed overlap
            n = len(gene_set)  # Sample size
            K = len(pathway_set)  # Population successes
            N = total_genes  # Population size
            
            # Calculate p-value
            p_value = stats.hypergeom.sf(k - 1, N, K, n)
            
            # Calculate expected overlap and fold enrichment
            expected = (n * K) / N
            fold_enrichment = k / expected if expected > 0 else float('inf')
            
            result = EnrichmentResult(
                term=pathway,
                observed=k,
                expected=expected,
                fold_enrichment=fold_enrichment,
                p_value=p_value,
                genes=list(overlap)
            )
            results.append(result)
        
        # Apply multiple testing correction
        if results:
            p_values = [r.p_value for r in results]
            _, corrected_p = self.corrector.correct_pvalues(p_values)
            
            for result, adj_p in zip(results, corrected_p):
                result.adjusted_p_value = adj_p
        
        # Sort by adjusted p-value
        results.sort(key=lambda x: x.adjusted_p_value or x.p_value)
        
        self.logger.info(f"Completed enrichment analysis for {len(results)} pathways")
        return results
    
    def binomial_enrichment(
        self,
        observed: List[int],
        expected: List[float],
        pathway_names: List[str]
    ) -> List[EnrichmentResult]:
        """
        Perform binomial enrichment test.
        
        Args:
            observed: Observed counts for each pathway
            expected: Expected counts for each pathway
            pathway_names: Names of pathways
            
        Returns:
            List of enrichment results
        """
        results = []
        
        for obs, exp, name in zip(observed, expected, pathway_names):
            if exp <= 0:
                continue
            
            # Binomial test
            p_value = stats.binom_test(obs, obs + exp, 0.5, alternative='greater')
            fold_enrichment = obs / exp if exp > 0 else float('inf')
            
            result = EnrichmentResult(
                term=name,
                observed=obs,
                expected=exp,
                fold_enrichment=fold_enrichment,
                p_value=p_value
            )
            results.append(result)
        
        # Apply multiple testing correction
        if results:
            p_values = [r.p_value for r in results]
            _, corrected_p = self.corrector.correct_pvalues(p_values)
            
            for result, adj_p in zip(results, corrected_p):
                result.adjusted_p_value = adj_p
        
        results.sort(key=lambda x: x.adjusted_p_value or x.p_value)
        return results


class StatisticalTester:
    """Comprehensive statistical testing framework."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def t_test(
        self, 
        group1: List[float], 
        group2: List[float],
        equal_var: bool = True,
        paired: bool = False
    ) -> StatisticalResult:
        """Perform t-test between two groups."""
        if paired:
            statistic, p_value = stats.ttest_rel(group1, group2)
            df = len(group1) - 1
        else:
            statistic, p_value = stats.ttest_ind(group1, group2, equal_var=equal_var)
            df = len(group1) + len(group2) - 2
        
        # Calculate effect size (Cohen's d)
        pooled_std = np.sqrt(((len(group1) - 1) * np.var(group1, ddof=1) + 
                             (len(group2) - 1) * np.var(group2, ddof=1)) / 
                            (len(group1) + len(group2) - 2))
        effect_size = (np.mean(group1) - np.mean(group2)) / pooled_std
        
        return StatisticalResult(
            test_name="t-test",
            statistic=statistic,
            p_value=p_value,
            effect_size=effect_size,
            degrees_freedom=df,
            sample_size=len(group1) + len(group2)
        )
    
    def mann_whitney_u(
        self, 
        group1: List[float], 
        group2: List[float],
        alternative: str = 'two-sided'
    ) -> StatisticalResult:
        """Perform Mann-Whitney U test."""
        statistic, p_value = mannwhitneyu(group1, group2, alternative=alternative)
        
        # Calculate effect size (rank-biserial correlation)
        n1, n2 = len(group1), len(group2)
        effect_size = 1 - (2 * statistic) / (n1 * n2)
        
        return StatisticalResult(
            test_name="Mann-Whitney U",
            statistic=statistic,
            p_value=p_value,
            effect_size=effect_size,
            sample_size=n1 + n2
        )
    
    def chi_square_test(
        self, 
        contingency_table: np.ndarray
    ) -> StatisticalResult:
        """Perform chi-square test of independence."""
        chi2, p_value, dof, expected = chi2_contingency(contingency_table)
        
        # Calculate Cramer's V (effect size)
        n = contingency_table.sum()
        min_dim = min(contingency_table.shape) - 1
        cramers_v = np.sqrt(chi2 / (n * min_dim))
        
        return StatisticalResult(
            test_name="Chi-square",
            statistic=chi2,
            p_value=p_value,
            effect_size=cramers_v,
            degrees_freedom=dof,
            sample_size=int(n)
        )
    
    def fisher_exact_test(
        self, 
        contingency_table: np.ndarray,
        alternative: str = 'two-sided'
    ) -> StatisticalResult:
        """Perform Fisher's exact test."""
        if contingency_table.shape != (2, 2):
            raise ValidationError("Fisher's exact test requires 2x2 contingency table")
        
        odds_ratio, p_value = fisher_exact(contingency_table, alternative=alternative)
        
        return StatisticalResult(
            test_name="Fisher's exact",
            statistic=odds_ratio,
            p_value=p_value,
            effect_size=np.log(odds_ratio) if odds_ratio > 0 else None,
            sample_size=int(contingency_table.sum())
        )
    
    def binomial_test(
        self, 
        successes: int, 
        trials: int, 
        expected_prob: float = 0.5,
        alternative: str = 'two-sided'
    ) -> StatisticalResult:
        """Perform binomial test."""
        p_value = stats.binom_test(successes, trials, expected_prob, alternative=alternative)
        
        # Effect size as difference from expected
        observed_prob = successes / trials
        effect_size = observed_prob - expected_prob
        
        return StatisticalResult(
            test_name="Binomial test",
            statistic=successes,
            p_value=p_value,
            effect_size=effect_size,
            sample_size=trials
        )


class CorrelationAnalyzer:
    """Analyze correlations and associations between variables."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def correlation_matrix(
        self, 
        data: pd.DataFrame, 
        method: str = 'pearson',
        min_periods: int = 1
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Calculate correlation matrix with p-values.
        
        Args:
            data: DataFrame with numeric columns
            method: Correlation method ('pearson', 'spearman', 'kendall')
            min_periods: Minimum number of observations for correlation
            
        Returns:
            Tuple of (correlation matrix, p-value matrix)
        """
        n_vars = len(data.columns)
        corr_matrix = np.zeros((n_vars, n_vars))
        p_matrix = np.zeros((n_vars, n_vars))
        
        for i, col1 in enumerate(data.columns):
            for j, col2 in enumerate(data.columns):
                if i == j:
                    corr_matrix[i, j] = 1.0
                    p_matrix[i, j] = 0.0
                else:
                    x = data[col1].dropna()
                    y = data[col2].dropna()
                    
                    # Find common indices
                    common_idx = x.index.intersection(y.index)
                    if len(common_idx) >= min_periods:
                        x_common = x.loc[common_idx]
                        y_common = y.loc[common_idx]
                        
                        if method == 'pearson':
                            corr, p_val = stats.pearsonr(x_common, y_common)
                        elif method == 'spearman':
                            corr, p_val = stats.spearmanr(x_common, y_common)
                        elif method == 'kendall':
                            corr, p_val = stats.kendalltau(x_common, y_common)
                        else:
                            raise ValidationError(f"Unknown correlation method: {method}")
                        
                        corr_matrix[i, j] = corr
                        p_matrix[i, j] = p_val
                    else:
                        corr_matrix[i, j] = np.nan
                        p_matrix[i, j] = np.nan
        
        corr_df = pd.DataFrame(corr_matrix, index=data.columns, columns=data.columns)
        p_df = pd.DataFrame(p_matrix, index=data.columns, columns=data.columns)
        
        return corr_df, p_df
    
    def partial_correlation(
        self, 
        data: pd.DataFrame, 
        x: str, 
        y: str, 
        control_vars: List[str]
    ) -> StatisticalResult:
        """
        Calculate partial correlation controlling for other variables.
        
        Args:
            data: DataFrame with variables
            x: First variable
            y: Second variable  
            control_vars: Variables to control for
            
        Returns:
            Statistical result with partial correlation
        """
        # Prepare data
        all_vars = [x, y] + control_vars
        clean_data = data[all_vars].dropna()
        
        if len(clean_data) < len(all_vars) + 2:
            raise ValidationError("Insufficient data for partial correlation")
        
        # Calculate correlation matrix
        corr_matrix = clean_data.corr().values
        
        # Calculate partial correlation using matrix inverse
        try:
            inv_corr = np.linalg.inv(corr_matrix)
            partial_corr = -inv_corr[0, 1] / np.sqrt(inv_corr[0, 0] * inv_corr[1, 1])
            
            # Calculate t-statistic and p-value
            n = len(clean_data)
            k = len(control_vars)
            df = n - k - 2
            
            t_stat = partial_corr * np.sqrt(df / (1 - partial_corr**2))
            p_value = 2 * (1 - stats.t.cdf(np.abs(t_stat), df))
            
            return StatisticalResult(
                test_name="Partial correlation",
                statistic=partial_corr,
                p_value=p_value,
                degrees_freedom=df,
                sample_size=n
            )
            
        except np.linalg.LinAlgError:
            raise ProcessingError("Cannot compute partial correlation: singular matrix")


class PowerAnalysis:
    """Power analysis for sample size determination and effect detection."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def power_t_test(
        self, 
        effect_size: float, 
        alpha: float = 0.05, 
        power: float = 0.8,
        alternative: str = 'two-sided'
    ) -> Dict[str, float]:
        """
        Calculate sample size for t-test given effect size and power.
        
        Args:
            effect_size: Cohen's d
            alpha: Type I error rate
            power: Desired statistical power
            alternative: Test type ('two-sided', 'one-sided')
            
        Returns:
            Dictionary with power analysis results
        """
        from statsmodels.stats.power import ttest_power
        
        # Calculate required sample size
        sample_size = ttest_power(
            effect_size=effect_size,
            power=power,
            alpha=alpha,
            alternative=alternative
        )
        
        return {
            'required_sample_size': int(np.ceil(sample_size)),
            'effect_size': effect_size,
            'alpha': alpha,
            'power': power,
            'alternative': alternative
        }
    
    def achieved_power(
        self, 
        effect_size: float, 
        sample_size: int, 
        alpha: float = 0.05,
        alternative: str = 'two-sided'
    ) -> float:
        """Calculate achieved power given effect size and sample size."""
        from statsmodels.stats.power import ttest_power
        
        return ttest_power(
            effect_size=effect_size,
            nobs=sample_size,
            alpha=alpha,
            alternative=alternative
        )


class StatisticalSummary:
    """Generate comprehensive statistical summaries."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def descriptive_stats(self, data: pd.DataFrame) -> pd.DataFrame:
        """Generate descriptive statistics for all numeric columns."""
        numeric_cols = data.select_dtypes(include=[np.number]).columns
        
        stats_dict = {}
        for col in numeric_cols:
            series = data[col].dropna()
            
            stats_dict[col] = {
                'count': len(series),
                'mean': series.mean(),
                'std': series.std(),
                'min': series.min(),
                'q25': series.quantile(0.25),
                'median': series.median(),
                'q75': series.quantile(0.75),
                'max': series.max(),
                'skewness': stats.skew(series),
                'kurtosis': stats.kurtosis(series),
                'normality_p': stats.shapiro(series)[1] if len(series) <= 5000 else stats.jarque_bera(series)[1]
            }
        
        return pd.DataFrame(stats_dict).T
    
    def group_comparison_summary(
        self, 
        data: pd.DataFrame, 
        group_col: str, 
        value_cols: List[str]
    ) -> pd.DataFrame:
        """Generate statistical comparison summary between groups."""
        tester = StatisticalTester()
        corrector = MultipleTestingCorrector()
        
        results = []
        groups = data[group_col].unique()
        
        if len(groups) != 2:
            raise ValidationError("Group comparison requires exactly 2 groups")
        
        group1_data = data[data[group_col] == groups[0]]
        group2_data = data[data[group_col] == groups[1]]
        
        for col in value_cols:
            if col not in data.columns:
                continue
            
            g1_values = group1_data[col].dropna().tolist()
            g2_values = group2_data[col].dropna().tolist()
            
            if len(g1_values) < 3 or len(g2_values) < 3:
                continue
            
            # Perform t-test and Mann-Whitney U test
            t_result = tester.t_test(g1_values, g2_values)
            mw_result = tester.mann_whitney_u(g1_values, g2_values)
            
            results.append({
                'variable': col,
                'group1_mean': np.mean(g1_values),
                'group1_std': np.std(g1_values),
                'group1_n': len(g1_values),
                'group2_mean': np.mean(g2_values),
                'group2_std': np.std(g2_values),
                'group2_n': len(g2_values),
                't_statistic': t_result.statistic,
                't_pvalue': t_result.p_value,
                'cohens_d': t_result.effect_size,
                'mw_statistic': mw_result.statistic,
                'mw_pvalue': mw_result.p_value,
                'rank_biserial_r': mw_result.effect_size
            })
        
        result_df = pd.DataFrame(results)
        
        # Apply multiple testing correction
        if not result_df.empty:
            _, t_corrected = corrector.correct_pvalues(result_df['t_pvalue'].tolist())
            _, mw_corrected = corrector.correct_pvalues(result_df['mw_pvalue'].tolist())
            
            result_df['t_pvalue_corrected'] = t_corrected
            result_df['mw_pvalue_corrected'] = mw_corrected
        
        return result_df