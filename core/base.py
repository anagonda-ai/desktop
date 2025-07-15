"""
Base classes for bioinformatics operations.

This module provides abstract base classes and mixins that define common interfaces
and functionality for all bioinformatics processors and components.
"""

import time
import logging
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Type, Union
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import threading
from dataclasses import dataclass

from .types import (
    AnalysisResult, 
    ProcessingStatus, 
    AnalysisType,
    ProcessingConfig,
    FilePath
)
from .exceptions import (
    BioinformaticsError, 
    ProcessingError, 
    ValidationError,
    ErrorCode
)
from .config import get_config


class LoggerMixin:
    """Mixin to provide logging capabilities to classes."""
    
    @property
    def logger(self) -> logging.Logger:
        """Get logger instance for this class."""
        if not hasattr(self, '_logger'):
            self._logger = logging.getLogger(self.__class__.__name__)
        return self._logger


class ValidatorMixin:
    """Mixin to provide validation capabilities."""
    
    def validate_file_exists(self, file_path: FilePath, file_type: str = "file") -> Path:
        """Validate that a file exists and is readable."""
        path = Path(file_path)
        if not path.exists():
            raise ValidationError(
                f"{file_type.capitalize()} not found: {path}",
                error_code=ErrorCode.FILE_NOT_FOUND,
                context={'file_path': str(path), 'file_type': file_type}
            )
        if not path.is_file():
            raise ValidationError(
                f"Path is not a file: {path}",
                error_code=ErrorCode.FILE_FORMAT_INVALID,
                context={'file_path': str(path)}
            )
        return path
    
    def validate_directory_exists(self, dir_path: FilePath, create: bool = False) -> Path:
        """Validate that a directory exists, optionally create it."""
        path = Path(dir_path)
        if not path.exists():
            if create:
                path.mkdir(parents=True, exist_ok=True)
                self.logger.info(f"Created directory: {path}")
            else:
                raise ValidationError(
                    f"Directory not found: {path}",
                    error_code=ErrorCode.FILE_NOT_FOUND,
                    context={'directory_path': str(path)}
                )
        return path
    
    def validate_numeric_range(
        self, 
        value: Union[int, float], 
        min_val: Optional[Union[int, float]] = None,
        max_val: Optional[Union[int, float]] = None,
        field_name: str = "value"
    ) -> Union[int, float]:
        """Validate that a numeric value is within specified range."""
        if min_val is not None and value < min_val:
            raise ValidationError(
                f"{field_name} must be >= {min_val}, got {value}",
                error_code=ErrorCode.INVALID_PARAMETERS,
                field_name=field_name,
                invalid_value=value
            )
        if max_val is not None and value > max_val:
            raise ValidationError(
                f"{field_name} must be <= {max_val}, got {value}",
                error_code=ErrorCode.INVALID_PARAMETERS,
                field_name=field_name,
                invalid_value=value
            )
        return value


class PerformanceMonitorMixin:
    """Mixin to provide performance monitoring capabilities."""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._performance_data: Dict[str, List[float]] = {}
        self._lock = threading.Lock()
    
    def _record_performance(self, operation: str, duration: float):
        """Record performance data for an operation."""
        with self._lock:
            if operation not in self._performance_data:
                self._performance_data[operation] = []
            self._performance_data[operation].append(duration)
    
    def get_performance_stats(self, operation: Optional[str] = None) -> Dict[str, Dict[str, float]]:
        """Get performance statistics."""
        with self._lock:
            if operation:
                data = self._performance_data.get(operation, [])
                if not data:
                    return {}
                return {
                    operation: {
                        'count': len(data),
                        'total': sum(data),
                        'average': sum(data) / len(data),
                        'min': min(data),
                        'max': max(data)
                    }
                }
            else:
                stats = {}
                for op, data in self._performance_data.items():
                    if data:
                        stats[op] = {
                            'count': len(data),
                            'total': sum(data),
                            'average': sum(data) / len(data),
                            'min': min(data),
                            'max': max(data)
                        }
                return stats


def performance_monitor(func):
    """Decorator to monitor function performance."""
    def wrapper(self, *args, **kwargs):
        start_time = time.time()
        try:
            result = func(self, *args, **kwargs)
            return result
        finally:
            duration = time.time() - start_time
            if hasattr(self, '_record_performance'):
                self._record_performance(func.__name__, duration)
            if hasattr(self, 'logger'):
                self.logger.debug(f"{func.__name__} completed in {duration:.3f}s")
    return wrapper


class BaseProcessor(LoggerMixin, ValidatorMixin, PerformanceMonitorMixin, ABC):
    """
    Abstract base class for all bioinformatics processors.
    
    Provides common functionality for validation, logging, configuration,
    and error handling that all processors should implement.
    """
    
    def __init__(self, config: Optional[ProcessingConfig] = None):
        super().__init__()
        self.config = config or ProcessingConfig()
        self.global_config = get_config()
        self._status = ProcessingStatus.PENDING
        self._results: Dict[str, Any] = {}
        self._errors: List[str] = []
    
    @property
    def status(self) -> ProcessingStatus:
        """Current processing status."""
        return self._status
    
    @property
    def results(self) -> Dict[str, Any]:
        """Processing results."""
        return self._results.copy()
    
    @property
    def errors(self) -> List[str]:
        """Processing errors."""
        return self._errors.copy()
    
    @abstractmethod
    def validate_input(self, data: Any) -> None:
        """
        Validate input data before processing.
        
        Args:
            data: Input data to validate
            
        Raises:
            ValidationError: If input data is invalid
        """
        pass
    
    @abstractmethod
    def process(self, data: Any, **kwargs) -> Any:
        """
        Process the input data.
        
        Args:
            data: Input data to process
            **kwargs: Additional processing parameters
            
        Returns:
            Processed results
            
        Raises:
            ProcessingError: If processing fails
        """
        pass
    
    @performance_monitor
    def run(self, data: Any, **kwargs) -> AnalysisResult:
        """
        Execute the complete processing pipeline.
        
        Args:
            data: Input data to process
            **kwargs: Additional processing parameters
            
        Returns:
            AnalysisResult with processing results and metadata
        """
        start_time = time.time()
        self._status = ProcessingStatus.IN_PROGRESS
        self._results.clear()
        self._errors.clear()
        
        try:
            self.logger.info(f"Starting {self.__class__.__name__} processing")
            
            # Validate input
            self.validate_input(data)
            
            # Process data
            results = self.process(data, **kwargs)
            
            # Store results
            self._results.update(results if isinstance(results, dict) else {'results': results})
            self._status = ProcessingStatus.COMPLETED
            
            processing_time = time.time() - start_time
            self.logger.info(f"Processing completed in {processing_time:.3f}s")
            
            return AnalysisResult(
                analysis_type=getattr(self, 'analysis_type', AnalysisType.BLAST_SEARCH),
                status=self._status,
                results=self._results,
                processing_time=processing_time,
                metadata={
                    'processor': self.__class__.__name__,
                    'config': self.config.__dict__,
                    'performance': self.get_performance_stats()
                }
            )
            
        except Exception as e:
            self._status = ProcessingStatus.FAILED
            error_msg = f"Processing failed: {str(e)}"
            self._errors.append(error_msg)
            self.logger.error(error_msg, exc_info=True)
            
            processing_time = time.time() - start_time
            
            return AnalysisResult(
                analysis_type=getattr(self, 'analysis_type', AnalysisType.BLAST_SEARCH),
                status=self._status,
                results=self._results,
                processing_time=processing_time,
                error_message=error_msg,
                metadata={
                    'processor': self.__class__.__name__,
                    'errors': self._errors
                }
            )


class BatchAPIProcessor(BaseProcessor):
    """
    Base class for processors that make batch API calls.
    
    Provides functionality for rate limiting, batch processing,
    and concurrent API operations.
    """
    
    def __init__(
        self, 
        config: Optional[ProcessingConfig] = None,
        batch_size: int = 10,
        rate_limit: float = 0.2,
        max_workers: int = 4
    ):
        super().__init__(config)
        self.batch_size = batch_size
        self.rate_limit = rate_limit
        self.max_workers = max_workers
        self._last_request_time = 0.0
        self._request_lock = threading.Lock()
    
    def _rate_limit_wait(self):
        """Enforce rate limiting between requests."""
        with self._request_lock:
            current_time = time.time()
            time_since_last = current_time - self._last_request_time
            if time_since_last < self.rate_limit:
                sleep_time = self.rate_limit - time_since_last
                time.sleep(sleep_time)
            self._last_request_time = time.time()
    
    def create_batches(self, items: List[Any]) -> List[List[Any]]:
        """Split items into batches for processing."""
        batches = []
        for i in range(0, len(items), self.batch_size):
            batch = items[i:i + self.batch_size]
            batches.append(batch)
        return batches
    
    @abstractmethod
    def process_batch(self, batch: List[Any], **kwargs) -> List[Any]:
        """
        Process a single batch of items.
        
        Args:
            batch: List of items to process
            **kwargs: Additional parameters
            
        Returns:
            List of processed results
        """
        pass
    
    @performance_monitor
    def process_batches_concurrent(
        self, 
        batches: List[List[Any]], 
        **kwargs
    ) -> List[Any]:
        """
        Process batches concurrently with rate limiting.
        
        Args:
            batches: List of batches to process
            **kwargs: Additional parameters
            
        Returns:
            Combined results from all batches
        """
        all_results = []
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all batches
            future_to_batch = {
                executor.submit(self._process_batch_with_rate_limit, batch, **kwargs): batch
                for batch in batches
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_batch):
                batch = future_to_batch[future]
                try:
                    batch_results = future.result()
                    all_results.extend(batch_results)
                    self.logger.debug(f"Completed batch with {len(batch)} items")
                except Exception as e:
                    self.logger.error(f"Batch processing failed: {e}")
                    raise ProcessingError(
                        f"Batch processing failed: {e}",
                        processor_name=self.__class__.__name__,
                        stage="batch_processing"
                    )
        
        return all_results
    
    def _process_batch_with_rate_limit(self, batch: List[Any], **kwargs) -> List[Any]:
        """Process a batch with rate limiting."""
        self._rate_limit_wait()
        return self.process_batch(batch, **kwargs)


class FileProcessor(BaseProcessor):
    """
    Base class for processors that handle file operations.
    
    Provides functionality for file validation, reading, writing,
    and format detection.
    """
    
    def __init__(self, config: Optional[ProcessingConfig] = None):
        super().__init__(config)
        self.supported_formats: List[str] = []
    
    def detect_file_format(self, file_path: FilePath) -> str:
        """
        Detect file format based on extension and content.
        
        Args:
            file_path: Path to the file
            
        Returns:
            Detected file format
            
        Raises:
            ValidationError: If format cannot be detected
        """
        path = Path(file_path)
        extension = path.suffix.lower()
        
        format_map = {
            '.fasta': 'fasta',
            '.fa': 'fasta',
            '.fas': 'fasta',
            '.csv': 'csv',
            '.tsv': 'tsv',
            '.txt': 'text',
            '.json': 'json',
            '.xml': 'xml',
            '.gff': 'gff',
            '.gff3': 'gff3',
            '.gtf': 'gtf',
            '.gb': 'genbank',
            '.gbk': 'genbank'
        }
        
        detected_format = format_map.get(extension)
        if not detected_format:
            raise ValidationError(
                f"Unable to detect file format for: {path}",
                error_code=ErrorCode.FILE_FORMAT_INVALID,
                context={'file_path': str(path), 'extension': extension}
            )
        
        return detected_format
    
    def validate_file_format(self, file_path: FilePath, expected_format: Optional[str] = None) -> str:
        """
        Validate file format against expected format.
        
        Args:
            file_path: Path to the file
            expected_format: Expected file format (if None, auto-detect)
            
        Returns:
            Validated file format
            
        Raises:
            ValidationError: If format validation fails
        """
        detected_format = self.detect_file_format(file_path)
        
        if expected_format and detected_format != expected_format:
            raise ValidationError(
                f"File format mismatch. Expected {expected_format}, got {detected_format}",
                error_code=ErrorCode.FILE_FORMAT_INVALID,
                context={
                    'file_path': str(file_path),
                    'expected': expected_format,
                    'detected': detected_format
                }
            )
        
        if self.supported_formats and detected_format not in self.supported_formats:
            raise ValidationError(
                f"Unsupported file format: {detected_format}. Supported: {self.supported_formats}",
                error_code=ErrorCode.FILE_FORMAT_INVALID,
                context={
                    'file_path': str(file_path),
                    'format': detected_format,
                    'supported': self.supported_formats
                }
            )
        
        return detected_format