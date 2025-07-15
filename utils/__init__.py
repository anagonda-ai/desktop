"""
Utility modules for bioinformatics operations.

This package provides shared utilities for file operations, API clients,
concurrent processing, logging, and visualization.
"""

from .file_operations import (
    CSVHandler,
    FASTAHandler,
    FileManager,
    SequenceFileReader,
    SequenceFileWriter,
    FileFormatDetector
)
from .api_clients import (
    KEGGClient,
    NCBIClient,
    UniProtClient,
    DatabaseClientFactory,
    APIResponse,
    CacheManager
)
from .concurrent import (
    BatchProcessor,
    RateLimitedExecutor,
    ProgressTracker,
    MemoryAwareBatchProcessor,
    ExecutorType,
    TaskResult
)
from .statistical import (
    MultipleTestingCorrector,
    EnrichmentAnalyzer,
    StatisticalTester,
    CorrelationAnalyzer,
    PowerAnalysis,
    StatisticalSummary,
    CorrectionMethod,
    StatisticalTest,
    StatisticalResult,
    EnrichmentResult
)
from .visualization import (
    ScientificPlotter,
    BioinformaticsThemes,
    PlotConfig,
    PlotTemplates,
    save_plot,
    PlotType,
    ColorScheme
)
from .logging import (
    setup_logging,
    get_logger,
    log_performance,
    log_workflow_start,
    log_workflow_end,
    BioinformaticsLogger
)

__all__ = [
    # File operations
    "CSVHandler",
    "FASTAHandler", 
    "FileManager",
    "SequenceFileReader",
    "SequenceFileWriter",
    "FileFormatDetector",
    # API clients
    "KEGGClient",
    "NCBIClient",
    "UniProtClient", 
    "DatabaseClientFactory",
    "APIResponse",
    "CacheManager",
    # Concurrent processing
    "BatchProcessor",
    "RateLimitedExecutor",
    "ProgressTracker",
    "MemoryAwareBatchProcessor",
    "ExecutorType",
    "TaskResult",
    # Statistical analysis
    "MultipleTestingCorrector",
    "EnrichmentAnalyzer",
    "StatisticalTester",
    "CorrelationAnalyzer",
    "PowerAnalysis",
    "StatisticalSummary",
    "CorrectionMethod",
    "StatisticalTest",
    "StatisticalResult",
    "EnrichmentResult",
    # Visualization
    "ScientificPlotter",
    "BioinformaticsThemes",
    "PlotConfig",
    "PlotTemplates",
    "save_plot",
    "PlotType",
    "ColorScheme",
    # Logging
    "setup_logging",
    "get_logger",
    "log_performance",
    "log_workflow_start",
    "log_workflow_end",
    "BioinformaticsLogger"
]