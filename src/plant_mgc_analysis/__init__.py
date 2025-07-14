"""
Plant MGC Analysis Pipeline
==========================

A comprehensive bioinformatics toolkit for identifying and analyzing 
biosynthetic gene clusters in plant genomes.

This package provides tools for:
- Metabolic gene cluster detection and analysis
- Genomic sliding window analysis
- BLAST-based homology searches
- Machine learning predictions for MGC classification
- Phylogenetic analysis of biosynthetic pathways
- Integration with multiple biological databases (KEGG, MiBIG, PlantCyc)

Modules:
    core: Core utilities and data structures
    genomics: Genomic analysis tools
    metabolic: Metabolic pathway analysis
    machine_learning: ML models for MGC prediction
    cli: Command-line interface
    config: Configuration management
    utils: Utility functions
    data: Data handling and processing

Example:
    >>> from plant_mgc_analysis import MGCAnalyzer
    >>> analyzer = MGCAnalyzer()
    >>> results = analyzer.analyze_genome("path/to/genome.fasta")
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("plant-mgc-analysis")
except PackageNotFoundError:
    __version__ = "unknown"

__author__ = "Itay Mayrose Lab"
__email__ = "itay.mayrose@example.com"
__license__ = "MIT"

# Core imports
from .core.analyzer import MGCAnalyzer
from .core.exceptions import (
    PlantMGCError,
    AnalysisError,
    ConfigurationError,
    DataError,
)
from .config.settings import get_settings

# Main API
__all__ = [
    "__version__",
    "MGCAnalyzer",
    "PlantMGCError",
    "AnalysisError", 
    "ConfigurationError",
    "DataError",
    "get_settings",
]