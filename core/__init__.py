"""
Core bioinformatics toolkit module.

This module provides the foundational classes and utilities for the bioinformatics
analysis pipeline, including base classes, configuration management, and common types.
"""

from .base import BaseProcessor, BatchAPIProcessor, FileProcessor
from .config import Config, DatabaseConfig, ComputeConfig
from .exceptions import (
    BioinformaticsError,
    ValidationError,
    FileOperationError,
    APIError,
    ProcessingError
)
from .types import (
    GenomeData,
    SequenceData,
    AnalysisResult,
    ProcessingStatus,
    FileFormat
)

__version__ = "1.0.0"
__author__ = "Itay Mayrose Lab"

__all__ = [
    "BaseProcessor",
    "BatchAPIProcessor", 
    "FileProcessor",
    "Config",
    "DatabaseConfig",
    "ComputeConfig",
    "BioinformaticsError",
    "ValidationError",
    "FileOperationError", 
    "APIError",
    "ProcessingError",
    "GenomeData",
    "SequenceData",
    "AnalysisResult",
    "ProcessingStatus",
    "FileFormat"
]