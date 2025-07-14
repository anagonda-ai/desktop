"""
Base classes and interfaces for Plant MGC Analysis Pipeline.

This module provides the foundational abstract base classes and interfaces
that define the structure and contracts for all analysis components.
"""

import time
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Protocol, TypeVar, Generic
from pathlib import Path
from dataclasses import dataclass
from contextlib import contextmanager

from loguru import logger

from .types import (
    AnalysisResult, 
    AnalysisType, 
    GenomeData,
    MGCCandidate,
    PathLike,
)
from .exceptions import AnalysisError, ValidationError
from ..config.settings import get_settings
from ..utils.logging import LoggerMixin, performance_monitor


T = TypeVar('T')
R = TypeVar('R')


class Validator(Protocol):
    """Protocol for data validation."""
    
    def validate(self, data: Any) -> bool:
        """Validate data and return True if valid."""
        ...
    
    def get_errors(self) -> List[str]:
        """Get list of validation errors."""
        ...


class Cacheable(Protocol):
    """Protocol for cacheable operations."""
    
    def get_cache_key(self) -> str:
        """Generate cache key for this operation."""
        ...
    
    def should_cache(self) -> bool:
        """Determine if this operation should be cached."""
        ...


class Processor(Protocol[T, R]):
    """Protocol for data processors."""
    
    def process(self, data: T) -> R:
        """Process input data and return result."""
        ...


@dataclass
class ProcessingContext:
    """Context for processing operations."""
    
    session_id: str
    temp_dir: Path
    config: Dict[str, Any]
    progress_callback: Optional[callable] = None
    
    def update_progress(self, message: str, percentage: Optional[float] = None) -> None:
        """Update progress if callback is available."""
        if self.progress_callback:
            self.progress_callback(message, percentage)


class BioinformaticsProcessor(LoggerMixin, ABC):
    """
    Abstract base class for all bioinformatics processors.
    
    Provides common functionality including configuration management,
    logging, error handling, and progress tracking.
    """
    
    def __init__(
        self,
        config: Optional[Dict[str, Any]] = None,
        session_id: Optional[str] = None,
    ):
        """
        Initialize processor with configuration.
        
        Args:
            config: Optional configuration dictionary
            session_id: Optional session identifier
        """
        self.settings = get_settings()
        if config:
            self.settings.update_config(**config)
        
        self.session_id = session_id or f"{self.__class__.__name__}_{int(time.time())}"
        self.temp_dir = self.settings.cache_dir / self.session_id
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        
        self._cache: Dict[str, Any] = {}
        self._processing_context: Optional[ProcessingContext] = None
        
        self.logger.info(f"Initialized {self.__class__.__name__} with session {self.session_id}")
    
    @property
    def context(self) -> ProcessingContext:
        """Get current processing context."""
        if self._processing_context is None:
            self._processing_context = ProcessingContext(
                session_id=self.session_id,
                temp_dir=self.temp_dir,
                config=self.settings.model_dump(),
            )
        return self._processing_context
    
    @contextmanager
    def processing_session(self, progress_callback: Optional[callable] = None):
        """Context manager for processing session."""
        self._processing_context = ProcessingContext(
            session_id=self.session_id,
            temp_dir=self.temp_dir,
            config=self.settings.model_dump(),
            progress_callback=progress_callback,
        )
        
        try:
            self.logger.info(f"Starting processing session {self.session_id}")
            yield self.context
        finally:
            self.logger.info(f"Ended processing session {self.session_id}")
            self._processing_context = None
    
    def cache_result(self, key: str, value: Any) -> None:
        """Cache a result with given key."""
        self._cache[key] = value
        self.logger.debug(f"Cached result with key: {key}")
    
    def get_cached_result(self, key: str) -> Optional[Any]:
        """Get cached result by key."""
        result = self._cache.get(key)
        if result is not None:
            self.logger.debug(f"Retrieved cached result for key: {key}")
        return result
    
    def clear_cache(self) -> None:
        """Clear all cached results."""
        self._cache.clear()
        self.logger.debug("Cleared cache")
    
    @abstractmethod
    def validate_input(self, data: Any) -> None:
        """
        Validate input data.
        
        Args:
            data: Input data to validate
            
        Raises:
            ValidationError: If validation fails
        """
        pass
    
    @abstractmethod
    def process(self, data: Any, **kwargs) -> Any:
        """
        Process input data.
        
        Args:
            data: Input data to process
            **kwargs: Additional processing parameters
            
        Returns:
            Processed data
        """
        pass
    
    def cleanup(self) -> None:
        """Clean up temporary resources."""
        if self.temp_dir.exists():
            import shutil
            shutil.rmtree(self.temp_dir)
            self.logger.info(f"Cleaned up temporary directory: {self.temp_dir}")


class AnalysisEngine(BioinformaticsProcessor):
    """
    Abstract base class for analysis engines.
    
    Provides structure for implementing specific analysis algorithms
    while maintaining consistent interface and error handling.
    """
    
    def __init__(
        self,
        analysis_type: AnalysisType,
        config: Optional[Dict[str, Any]] = None,
        session_id: Optional[str] = None,
    ):
        """
        Initialize analysis engine.
        
        Args:
            analysis_type: Type of analysis this engine performs
            config: Optional configuration dictionary
            session_id: Optional session identifier
        """
        super().__init__(config, session_id)
        self.analysis_type = analysis_type
        self._results: List[AnalysisResult] = []
    
    @abstractmethod
    def analyze(
        self, 
        genome_data: GenomeData, 
        parameters: Optional[Dict[str, Any]] = None
    ) -> AnalysisResult:
        """
        Perform analysis on genome data.
        
        Args:
            genome_data: Genome data to analyze
            parameters: Optional analysis parameters
            
        Returns:
            Analysis result
        """
        pass
    
    @performance_monitor
    def run_analysis(
        self,
        genome_data: GenomeData,
        parameters: Optional[Dict[str, Any]] = None,
        progress_callback: Optional[callable] = None,
    ) -> AnalysisResult:
        """
        Run complete analysis with error handling and logging.
        
        Args:
            genome_data: Genome data to analyze
            parameters: Optional analysis parameters
            progress_callback: Optional progress callback
            
        Returns:
            Analysis result
            
        Raises:
            AnalysisError: If analysis fails
        """
        with self.processing_session(progress_callback) as context:
            context.update_progress(f"Starting {self.analysis_type.value} analysis", 0)
            
            try:
                # Validate inputs
                self.validate_input(genome_data)
                if parameters:
                    self.validate_parameters(parameters)
                
                context.update_progress("Running analysis", 25)
                
                # Run analysis
                result = self.analyze(genome_data, parameters)
                
                context.update_progress("Analysis completed", 100)
                
                # Store result
                self._results.append(result)
                
                self.logger.info(
                    f"Analysis {self.analysis_type.value} completed successfully"
                )
                
                return result
                
            except Exception as e:
                error_msg = f"Analysis {self.analysis_type.value} failed: {str(e)}"
                self.logger.error(error_msg)
                
                # Create failed result
                result = AnalysisResult(
                    analysis_id=f"{self.session_id}_{self.analysis_type.value}",
                    analysis_type=self.analysis_type,
                    input_data={"genome": str(genome_data.genome_file)},
                    results={},
                    parameters=parameters or {},
                    status="failed",
                    error_message=str(e),
                )
                
                self._results.append(result)
                raise AnalysisError(error_msg, analysis_type=self.analysis_type.value)
    
    def validate_parameters(self, parameters: Dict[str, Any]) -> None:
        """
        Validate analysis parameters.
        
        Args:
            parameters: Parameters to validate
            
        Raises:
            ValidationError: If parameters are invalid
        """
        # Base implementation - override in subclasses for specific validation
        if not isinstance(parameters, dict):
            raise ValidationError(
                "Parameters must be a dictionary",
                field_name="parameters",
                field_value=parameters,
                expected_type=dict,
            )
    
    def get_results(self) -> List[AnalysisResult]:
        """Get all analysis results."""
        return self._results.copy()
    
    def get_latest_result(self) -> Optional[AnalysisResult]:
        """Get the most recent analysis result."""
        return self._results[-1] if self._results else None


class DatabaseInterface(BioinformaticsProcessor):
    """
    Abstract base class for database interfaces.
    
    Provides common functionality for database connections,
    query execution, and result processing.
    """
    
    def __init__(
        self,
        database_name: str,
        connection_params: Dict[str, Any],
        config: Optional[Dict[str, Any]] = None,
        session_id: Optional[str] = None,
    ):
        """
        Initialize database interface.
        
        Args:
            database_name: Name of the database
            connection_params: Database connection parameters
            config: Optional configuration dictionary
            session_id: Optional session identifier
        """
        super().__init__(config, session_id)
        self.database_name = database_name
        self.connection_params = connection_params
        self._connection = None
    
    @abstractmethod
    def connect(self) -> None:
        """Establish database connection."""
        pass
    
    @abstractmethod
    def disconnect(self) -> None:
        """Close database connection."""
        pass
    
    @abstractmethod
    def query(self, query: str, parameters: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
        """
        Execute query and return results.
        
        Args:
            query: Query string
            parameters: Optional query parameters
            
        Returns:
            Query results
        """
        pass
    
    @contextmanager
    def connection(self):
        """Context manager for database connection."""
        try:
            self.connect()
            yield self._connection
        finally:
            self.disconnect()
    
    def validate_input(self, data: Any) -> None:
        """Validate database query input."""
        if not isinstance(data, (str, dict)):
            raise ValidationError(
                "Database input must be a query string or parameters dict",
                field_name="query_input",
                field_value=type(data).__name__,
            )


class FileProcessor(BioinformaticsProcessor):
    """
    Abstract base class for file processors.
    
    Provides common functionality for file I/O operations,
    format validation, and batch processing.
    """
    
    def __init__(
        self,
        supported_formats: List[str],
        config: Optional[Dict[str, Any]] = None,
        session_id: Optional[str] = None,
    ):
        """
        Initialize file processor.
        
        Args:
            supported_formats: List of supported file formats
            config: Optional configuration dictionary
            session_id: Optional session identifier
        """
        super().__init__(config, session_id)
        self.supported_formats = supported_formats
    
    def validate_input(self, file_path: PathLike) -> None:
        """
        Validate input file.
        
        Args:
            file_path: Path to file
            
        Raises:
            ValidationError: If file is invalid
        """
        path = Path(file_path)
        
        if not path.exists():
            raise ValidationError(
                f"File does not exist: {path}",
                field_name="file_path",
                field_value=str(path),
            )
        
        if not path.is_file():
            raise ValidationError(
                f"Path is not a file: {path}",
                field_name="file_path",
                field_value=str(path),
            )
        
        # Check file format
        file_extension = path.suffix.lower()
        if file_extension not in self.supported_formats:
            raise ValidationError(
                f"Unsupported file format: {file_extension}. "
                f"Supported formats: {self.supported_formats}",
                field_name="file_format",
                field_value=file_extension,
            )
    
    @abstractmethod
    def read_file(self, file_path: PathLike) -> Any:
        """
        Read file and return data.
        
        Args:
            file_path: Path to file
            
        Returns:
            File data
        """
        pass
    
    @abstractmethod
    def write_file(self, data: Any, file_path: PathLike) -> None:
        """
        Write data to file.
        
        Args:
            data: Data to write
            file_path: Output file path
        """
        pass
    
    def process_batch(
        self,
        file_paths: List[PathLike],
        output_dir: PathLike,
        parallel: bool = True,
        max_workers: Optional[int] = None,
    ) -> List[Any]:
        """
        Process multiple files in batch.
        
        Args:
            file_paths: List of input file paths
            output_dir: Output directory
            parallel: Whether to process in parallel
            max_workers: Maximum number of worker processes
            
        Returns:
            List of processing results
        """
        from concurrent.futures import ProcessPoolExecutor, as_completed
        
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        if not parallel or len(file_paths) == 1:
            # Sequential processing
            results = []
            for file_path in file_paths:
                result = self.process(file_path)
                results.append(result)
            return results
        
        # Parallel processing
        max_workers = max_workers or self.settings.compute.max_workers
        results = []
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_file = {
                executor.submit(self.process, file_path): file_path
                for file_path in file_paths
            }
            
            for future in as_completed(future_to_file):
                file_path = future_to_file[future]
                try:
                    result = future.result()
                    results.append(result)
                    self.logger.info(f"Processed file: {file_path}")
                except Exception as e:
                    self.logger.error(f"Error processing file {file_path}: {e}")
                    results.append(None)
        
        return results


class StatisticalAnalyzer(BioinformaticsProcessor):
    """
    Abstract base class for statistical analysis components.
    
    Provides common statistical operations and hypothesis testing
    functionality for bioinformatics analyses.
    """
    
    def __init__(
        self,
        significance_threshold: float = 0.05,
        multiple_testing_correction: str = "fdr_bh",
        config: Optional[Dict[str, Any]] = None,
        session_id: Optional[str] = None,
    ):
        """
        Initialize statistical analyzer.
        
        Args:
            significance_threshold: P-value threshold for significance
            multiple_testing_correction: Method for multiple testing correction
            config: Optional configuration dictionary
            session_id: Optional session identifier
        """
        super().__init__(config, session_id)
        self.significance_threshold = significance_threshold
        self.multiple_testing_correction = multiple_testing_correction
    
    @abstractmethod
    def calculate_statistics(self, data: Any) -> Dict[str, float]:
        """
        Calculate statistical measures for data.
        
        Args:
            data: Input data
            
        Returns:
            Dictionary of statistical measures
        """
        pass
    
    @abstractmethod
    def perform_hypothesis_test(
        self,
        sample1: Any,
        sample2: Any,
        test_type: str = "two_sided",
    ) -> Dict[str, float]:
        """
        Perform hypothesis test between two samples.
        
        Args:
            sample1: First sample
            sample2: Second sample
            test_type: Type of test (two_sided, greater, less)
            
        Returns:
            Test results including p-value and test statistic
        """
        pass
    
    def apply_multiple_testing_correction(self, p_values: List[float]) -> List[float]:
        """
        Apply multiple testing correction to p-values.
        
        Args:
            p_values: List of uncorrected p-values
            
        Returns:
            List of corrected p-values
        """
        from scipy.stats import false_discovery_control
        
        if self.multiple_testing_correction == "fdr_bh":
            return false_discovery_control(p_values).tolist()
        else:
            # Add other correction methods as needed
            self.logger.warning(
                f"Unknown correction method: {self.multiple_testing_correction}. "
                "Using FDR-BH."
            )
            return false_discovery_control(p_values).tolist()
    
    def is_significant(self, p_value: float) -> bool:
        """
        Check if p-value is significant.
        
        Args:
            p_value: P-value to check
            
        Returns:
            True if significant
        """
        return p_value < self.significance_threshold