"""
Advanced visualization framework for bioinformatics data.

This module provides standardized, publication-ready plotting capabilities
with scientific themes, consistent styling, and automated layout optimization.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
from typing import Any, Dict, List, Optional, Tuple, Union
from pathlib import Path
from dataclasses import dataclass, field
from enum import Enum
import colorcet as cc
from scipy import stats
import warnings

from ..core.types import FilePath
from ..core.exceptions import ValidationError


class PlotType(Enum):
    """Supported plot types."""
    SCATTER = "scatter"
    LINE = "line"
    BAR = "bar"
    HEATMAP = "heatmap"
    HISTOGRAM = "histogram"
    BOXPLOT = "boxplot"
    VIOLIN = "violin"
    TREE = "tree"
    GENOME_TRACK = "genome_track"
    ENRICHMENT = "enrichment"
    CORRELATION = "correlation"
    MANHATTAN = "manhattan"
    VOLCANO = "volcano"


class ColorScheme(Enum):
    """Scientific color schemes."""
    VIRIDIS = "viridis"
    PLASMA = "plasma"
    CIVIDIS = "cividis"
    COOLWARM = "coolwarm"
    SPECTRAL = "Spectral"
    CATEGORY10 = "Category10"
    SET1 = "Set1"
    PAIRED = "Paired"
    FIRE = "fire"
    BMGY = "bmgy"


@dataclass
class PlotConfig:
    """Configuration for plot appearance and behavior."""
    width: int = 10
    height: int = 6
    dpi: int = 300
    style: str = "whitegrid"
    color_scheme: ColorScheme = ColorScheme.VIRIDIS
    font_size: int = 12
    title_size: int = 14
    label_size: int = 11
    tick_size: int = 10
    legend_size: int = 10
    save_formats: List[str] = field(default_factory=lambda: ["png", "svg"])
    transparent: bool = False
    tight_layout: bool = True
    
    # Statistical annotation settings
    show_statistics: bool = True
    alpha: float = 0.05
    correction_method: str = "fdr_bh"
    
    # Interactive plot settings
    interactive: bool = False
    show_hover: bool = True
    
    # Publication settings
    publication_ready: bool = True
    remove_spines: bool = True
    grid_alpha: float = 0.3


class BioinformaticsThemes:
    """
    Scientific themes and color palettes for bioinformatics visualizations.
    """
    
    @staticmethod
    def get_color_palette(scheme: ColorScheme, n_colors: Optional[int] = None) -> List[str]:
        """Get color palette for specified scheme."""
        color_maps = {
            ColorScheme.VIRIDIS: cc.glasbey_category10,
            ColorScheme.PLASMA: cc.fire,
            ColorScheme.CIVIDIS: cc.bmgy,
            ColorScheme.COOLWARM: ["#3498db", "#e74c3c"],
            ColorScheme.SPECTRAL: cc.rainbow,
            ColorScheme.CATEGORY10: px.colors.qualitative.Plotly,
            ColorScheme.SET1: px.colors.qualitative.Set1,
            ColorScheme.PAIRED: px.colors.qualitative.Paired,
            ColorScheme.FIRE: cc.fire,
            ColorScheme.BMGY: cc.bmgy
        }
        
        colors = color_maps.get(scheme, cc.glasbey_category10)
        if n_colors:
            if len(colors) >= n_colors:
                return colors[:n_colors]
            else:
                # Extend palette by cycling
                extended = []
                for i in range(n_colors):
                    extended.append(colors[i % len(colors)])
                return extended
        return colors
    
    @staticmethod
    def apply_publication_theme():
        """Apply publication-ready theme to matplotlib."""
        plt.rcParams.update({
            'font.size': 12,
            'axes.titlesize': 14,
            'axes.labelsize': 12,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 11,
            'figure.titlesize': 16,
            'font.family': 'sans-serif',
            'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans'],
            'axes.spines.top': False,
            'axes.spines.right': False,
            'axes.grid': True,
            'grid.alpha': 0.3,
            'lines.linewidth': 2,
            'patch.linewidth': 0.5,
            'axes.linewidth': 1.2,
            'xtick.major.width': 1.2,
            'ytick.major.width': 1.2,
            'xtick.minor.width': 0.8,
            'ytick.minor.width': 0.8,
            'savefig.dpi': 300,
            'savefig.bbox': 'tight',
            'savefig.pad_inches': 0.1
        })


class ScientificPlotter:
    """
    Advanced plotting class for scientific visualizations.
    
    Provides publication-ready plots with statistical annotations,
    consistent styling, and automated layout optimization.
    """
    
    def __init__(self, config: Optional[PlotConfig] = None):
        self.config = config or PlotConfig()
        self._setup_style()
        self.themes = BioinformaticsThemes()
    
    def _setup_style(self):
        """Setup plotting style based on configuration."""
        if self.config.publication_ready:
            self.themes.apply_publication_theme()
        
        sns.set_style(self.config.style)
        plt.rcParams['figure.figsize'] = (self.config.width, self.config.height)
        plt.rcParams['savefig.dpi'] = self.config.dpi
    
    def create_scatter_plot(
        self,
        data: pd.DataFrame,
        x: str,
        y: str,
        color: Optional[str] = None,
        size: Optional[str] = None,
        title: Optional[str] = None,
        **kwargs
    ) -> plt.Figure:
        """
        Create a scatter plot with optional grouping and statistical annotations.
        
        Args:
            data: DataFrame containing the data
            x: Column name for x-axis
            y: Column name for y-axis
            color: Column name for color grouping
            size: Column name for size mapping
            title: Plot title
            **kwargs: Additional plotting parameters
            
        Returns:
            matplotlib Figure object
        """
        fig, ax = plt.subplots(figsize=(self.config.width, self.config.height))
        
        # Create scatter plot
        if color:
            groups = data[color].unique()
            colors = self.themes.get_color_palette(self.config.color_scheme, len(groups))
            for i, group in enumerate(groups):
                group_data = data[data[color] == group]
                ax.scatter(
                    group_data[x], 
                    group_data[y],
                    c=colors[i],
                    label=group,
                    alpha=0.7,
                    s=60 if not size else group_data[size] * 10
                )
            ax.legend()
        else:
            ax.scatter(data[x], data[y], alpha=0.7, s=60)
        
        # Add correlation statistics if requested
        if self.config.show_statistics and not color:
            corr, p_value = stats.pearsonr(data[x].dropna(), data[y].dropna())
            ax.text(
                0.05, 0.95, 
                f'r = {corr:.3f}, p = {p_value:.3e}',
                transform=ax.transAxes,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8)
            )
        
        # Formatting
        ax.set_xlabel(x.replace('_', ' ').title())
        ax.set_ylabel(y.replace('_', ' ').title())
        if title:
            ax.set_title(title)
        
        if self.config.tight_layout:
            plt.tight_layout()
        
        return fig
    
    def create_heatmap(
        self,
        data: Union[pd.DataFrame, np.ndarray],
        title: Optional[str] = None,
        row_labels: Optional[List[str]] = None,
        col_labels: Optional[List[str]] = None,
        annotate: bool = True,
        **kwargs
    ) -> plt.Figure:
        """
        Create a heatmap with dendrograms and statistical annotations.
        
        Args:
            data: Data matrix for heatmap
            title: Plot title
            row_labels: Labels for rows
            col_labels: Labels for columns
            annotate: Whether to annotate cells with values
            **kwargs: Additional seaborn heatmap parameters
            
        Returns:
            matplotlib Figure object
        """
        fig, ax = plt.subplots(figsize=(self.config.width, self.config.height))
        
        # Create heatmap
        sns.heatmap(
            data,
            ax=ax,
            cmap=self.config.color_scheme.value,
            annot=annotate,
            fmt='.2f' if annotate else '',
            xticklabels=col_labels or True,
            yticklabels=row_labels or True,
            cbar_kws={'shrink': 0.8},
            **kwargs
        )
        
        if title:
            ax.set_title(title)
        
        plt.tight_layout()
        return fig
    
    def create_enrichment_plot(
        self,
        data: pd.DataFrame,
        term_col: str,
        pvalue_col: str,
        fold_change_col: str,
        title: Optional[str] = None,
        max_terms: int = 20
    ) -> plt.Figure:
        """
        Create an enrichment plot showing pathway/GO term enrichment.
        
        Args:
            data: DataFrame with enrichment results
            term_col: Column name for terms/pathways
            pvalue_col: Column name for p-values
            fold_change_col: Column name for fold changes
            title: Plot title
            max_terms: Maximum number of terms to display
            
        Returns:
            matplotlib Figure object
        """
        # Prepare data
        plot_data = data.copy()
        plot_data['-log10(p)'] = -np.log10(plot_data[pvalue_col])
        plot_data = plot_data.nlargest(max_terms, '-log10(p)')
        
        fig, ax = plt.subplots(figsize=(self.config.width, self.config.height))
        
        # Create bubble plot
        scatter = ax.scatter(
            plot_data[fold_change_col],
            range(len(plot_data)),
            s=plot_data['-log10(p)'] * 20,
            c=plot_data['-log10(p)'],
            cmap=self.config.color_scheme.value,
            alpha=0.7,
            edgecolors='black',
            linewidth=0.5
        )
        
        # Add significance threshold line
        if self.config.show_statistics:
            significance_threshold = -np.log10(self.config.alpha)
            ax.axhline(
                y=significance_threshold, 
                color='red', 
                linestyle='--', 
                alpha=0.5,
                label=f'p = {self.config.alpha}'
            )
        
        # Formatting
        ax.set_yticks(range(len(plot_data)))
        ax.set_yticklabels(plot_data[term_col])
        ax.set_xlabel('Fold Change')
        ax.set_ylabel('Pathways/Terms')
        if title:
            ax.set_title(title)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('-log10(p-value)')
        
        plt.tight_layout()
        return fig
    
    def create_volcano_plot(
        self,
        data: pd.DataFrame,
        log2fc_col: str,
        pvalue_col: str,
        gene_col: Optional[str] = None,
        title: Optional[str] = None,
        fc_threshold: float = 1.0,
        label_top: int = 10
    ) -> plt.Figure:
        """
        Create a volcano plot for differential expression analysis.
        
        Args:
            data: DataFrame with DE results
            log2fc_col: Column name for log2 fold changes
            pvalue_col: Column name for p-values
            gene_col: Column name for gene names (for labeling)
            title: Plot title
            fc_threshold: Fold change threshold for significance
            label_top: Number of top genes to label
            
        Returns:
            matplotlib Figure object
        """
        # Prepare data
        plot_data = data.copy()
        plot_data['-log10(p)'] = -np.log10(plot_data[pvalue_col])
        
        # Define significance categories
        plot_data['significance'] = 'Not Significant'
        sig_mask = (
            (np.abs(plot_data[log2fc_col]) > fc_threshold) & 
            (plot_data[pvalue_col] < self.config.alpha)
        )
        plot_data.loc[sig_mask, 'significance'] = 'Significant'
        
        fig, ax = plt.subplots(figsize=(self.config.width, self.config.height))
        
        # Create scatter plot with color coding
        colors = ['lightgray', 'red']
        for i, (significance, group_data) in enumerate(plot_data.groupby('significance')):
            ax.scatter(
                group_data[log2fc_col],
                group_data['-log10(p)'],
                c=colors[i],
                alpha=0.6,
                label=significance,
                s=30
            )
        
        # Add threshold lines
        ax.axhline(
            y=-np.log10(self.config.alpha), 
            color='blue', 
            linestyle='--', 
            alpha=0.5,
            label=f'p = {self.config.alpha}'
        )
        ax.axvline(x=fc_threshold, color='blue', linestyle='--', alpha=0.5)
        ax.axvline(x=-fc_threshold, color='blue', linestyle='--', alpha=0.5)
        
        # Label top genes
        if gene_col and label_top > 0:
            top_genes = plot_data.nlargest(label_top, '-log10(p)')
            for _, gene in top_genes.iterrows():
                ax.annotate(
                    gene[gene_col],
                    (gene[log2fc_col], gene['-log10(p)']),
                    xytext=(5, 5),
                    textcoords='offset points',
                    fontsize=8,
                    alpha=0.8
                )
        
        # Formatting
        ax.set_xlabel('log2(Fold Change)')
        ax.set_ylabel('-log10(p-value)')
        if title:
            ax.set_title(title)
        ax.legend()
        
        plt.tight_layout()
        return fig
    
    def create_genome_track(
        self,
        genes: List[Dict[str, Any]],
        chromosome: str,
        start: int,
        end: int,
        title: Optional[str] = None
    ) -> plt.Figure:
        """
        Create a genome track visualization showing genes and features.
        
        Args:
            genes: List of gene dictionaries with start, end, strand, name
            chromosome: Chromosome name
            start: Start coordinate
            end: End coordinate
            title: Plot title
            
        Returns:
            matplotlib Figure object
        """
        fig, ax = plt.subplots(figsize=(self.config.width, 4))
        
        # Filter genes in the region
        region_genes = [
            gene for gene in genes 
            if gene['start'] >= start and gene['end'] <= end
        ]
        
        # Draw genes
        colors = self.themes.get_color_palette(self.config.color_scheme, 2)
        for i, gene in enumerate(region_genes):
            y_pos = 0.5 + (i % 2) * 0.3  # Alternate gene positions
            
            # Gene body
            rect = patches.Rectangle(
                (gene['start'], y_pos - 0.1),
                gene['end'] - gene['start'],
                0.2,
                linewidth=1,
                edgecolor='black',
                facecolor=colors[0] if gene.get('strand', '+') == '+' else colors[1],
                alpha=0.7
            )
            ax.add_patch(rect)
            
            # Gene name
            ax.text(
                (gene['start'] + gene['end']) / 2,
                y_pos + 0.15,
                gene.get('name', f"Gene_{i}"),
                ha='center',
                va='bottom',
                fontsize=8,
                rotation=45 if len(gene.get('name', '')) > 10 else 0
            )
            
            # Strand arrow
            arrow_x = gene['end'] if gene.get('strand', '+') == '+' else gene['start']
            arrow_dx = 100 if gene.get('strand', '+') == '+' else -100
            ax.arrow(
                arrow_x, y_pos,
                arrow_dx, 0,
                head_width=0.05,
                head_length=50,
                fc='black',
                ec='black',
                alpha=0.8
            )
        
        # Formatting
        ax.set_xlim(start, end)
        ax.set_ylim(0, 1.5)
        ax.set_xlabel(f'{chromosome} Position (bp)')
        ax.set_ylabel('')
        ax.set_yticks([])
        
        if title:
            ax.set_title(title)
        
        plt.tight_layout()
        return fig
    
    def save_plot(
        self, 
        fig: plt.Figure, 
        filename: FilePath, 
        formats: Optional[List[str]] = None
    ):
        """
        Save plot in multiple formats with consistent settings.
        
        Args:
            fig: matplotlib Figure to save
            filename: Base filename (without extension)
            formats: List of formats to save (defaults to config setting)
        """
        formats = formats or self.config.save_formats
        base_path = Path(filename)
        
        for fmt in formats:
            output_path = base_path.with_suffix(f'.{fmt}')
            fig.savefig(
                output_path,
                format=fmt,
                dpi=self.config.dpi,
                transparent=self.config.transparent,
                bbox_inches='tight' if self.config.tight_layout else None,
                pad_inches=0.1
            )


def save_plot(
    fig: plt.Figure, 
    filename: FilePath, 
    config: Optional[PlotConfig] = None,
    formats: Optional[List[str]] = None
):
    """
    Convenience function to save plots with standard settings.
    
    Args:
        fig: matplotlib Figure to save
        filename: Base filename (without extension)
        config: Plot configuration
        formats: List of formats to save
    """
    plotter = ScientificPlotter(config)
    plotter.save_plot(fig, filename, formats)


# Example usage and templates
class PlotTemplates:
    """Ready-to-use plot templates for common bioinformatics visualizations."""
    
    @staticmethod
    def mgc_enrichment_plot(enrichment_data: pd.DataFrame) -> plt.Figure:
        """Template for MGC pathway enrichment visualization."""
        plotter = ScientificPlotter()
        return plotter.create_enrichment_plot(
            enrichment_data,
            term_col='pathway',
            pvalue_col='pvalue',
            fold_change_col='fold_change',
            title='MGC Pathway Enrichment Analysis'
        )
    
    @staticmethod
    def sliding_window_results(window_data: pd.DataFrame) -> plt.Figure:
        """Template for sliding window analysis results."""
        plotter = ScientificPlotter()
        return plotter.create_scatter_plot(
            window_data,
            x='position',
            y='enrichment_score',
            color='significant',
            title='Sliding Window Enrichment Analysis'
        )
    
    @staticmethod
    def correlation_matrix(correlation_data: pd.DataFrame) -> plt.Figure:
        """Template for correlation matrix visualization."""
        plotter = ScientificPlotter()
        return plotter.create_heatmap(
            correlation_data,
            title='Gene Expression Correlation Matrix',
            annotate=False
        )