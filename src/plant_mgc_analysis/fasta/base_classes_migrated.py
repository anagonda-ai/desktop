"""
Industry-Level Base Classes for Python Scripts Bioinformatics Pipeline.

This module provides foundational OOP classes for all bioinformatics scripts
with comprehensive error handling, logging, validation, and performance optimization.
"""

from ..config.settings import get_settings
from ..core.base import BioinformaticsProcessor
from ..core.exceptions import ProcessingError, ValidationError
from ..utils.file_operations import file_manager
from loguru import logger
import os
import sys
import time
import json
import pickle
import logging
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Any, Union, Type, Tuple
from dataclasses import dataclass, field
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from contextlib import contextmanager
import tempfile
import shutil

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@dataclass
class ProcessingResult:
    """Result container for processing operations."""
    
    success: bool
    data: Any = None
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    processing_time: float = 0.0
    
    def add_error(self, error: str) -> None:
        """Add error message."""
        self.errors.append(error)
        self.success = False
    
    def add_warning(self, warning: str) -> None:
        """Add warning message."""
        self.warnings.append(warning)
    
    def merge(self, other: 'ProcessingResult') -> 'ProcessingResult':
        """Merge with another result."""
        return ProcessingResult(
            success=self.success and other.success,
            data=other.data if other.data is not None else self.data,
            errors=self.errors + other.errors,
            warnings=self.warnings + other.warnings,
            metadata={**self.metadata, **other.metadata},
            processing_time=self.processing_time + other.processing_time
        )


@dataclass
class BioinformaticsConfig:
    """Configuration for bioinformatics operations."""
    
    input_dir: Path
    output_dir: Path
    temp_dir: Optional[Path] = None
    max_workers: int = 4
    chunk_size: int = 1000
    timeout: int = 3600
    retry_attempts: int = 3
    log_level: str = "INFO"
    cache_enabled: bool = True
    
    def __post_init__(self):
        """Validate and setup configuration."""
        # Ensure directories exist
        self.input_dir = Path(self.input_dir)
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        if self.temp_dir is None:
            self.temp_dir = Path(tempfile.mkdtemp(prefix="bioinformatics_"))
        else:
            self.temp_dir = Path(self.temp_dir)
            self.temp_dir.mkdir(parents=True, exist_ok=True)


class BioinformaticsLogger:
    """Enhanced logging system for bioinformatics operations."""
    
    def __init__(self, name: str, log_level: str = "INFO"):
        """Initialize logger."""
        self.logger = logging.getLogger(name)
        self.logger.setLevel(getattr(logging, log_level.upper()))
        
        if not self.logger.handlers:
            # Console handler
            console_handler = logging.StreamHandler()
            console_formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            console_handler.setFormatter(console_formatter)
            self.logger.addHandler(console_handler)
    
    def info(self, message: str, **kwargs) -> None:
        """Log info message."""
        self.logger.info(f"{message} {self._format_kwargs(kwargs)}")
    
    def warning(self, message: str, **kwargs) -> None:
        """Log warning message."""
        self.logger.warning(f"{message} {self._format_kwargs(kwargs)}")
    
    def error(self, message: str, **kwargs) -> None:
        """Log error message."""
        self.logger.error(f"{message} {self._format_kwargs(kwargs)}")
    
    def debug(self, message: str, **kwargs) -> None:
        """Log debug message."""
        self.logger.debug(f"{message} {self._format_kwargs(kwargs)}")
    
    def _format_kwargs(self, kwargs: Dict[str, Any]) -> str:
        """Format kwargs for logging."""
        if not kwargs:
            return ""
        formatted = []
        for key, value in kwargs.items():
            formatted.append(f"{key}={value}")
        return f"[{', '.join(formatted)}]"


class BaseProcessor(ABC):
    """
    Abstract base class for all bioinformatics processors.
    
    Provides common functionality for file handling, logging, error management,
    and performance monitoring.
    """
    
    def __init__(self, config: BioinformaticsConfig, name: Optional[str] = None):
        """Initialize processor."""
        self.config = config
        self.name = name or self.__class__.__name__
        self.logger = BioinformaticsLogger(self.name, config.log_level)
        
        # Performance tracking
        self._start_time: Optional[float] = None
        self._stats = {
            "files_processed": 0,
            "errors_encountered": 0,
            "total_processing_time": 0.0,
        }
        
        # Cache
        self._cache: Dict[str, Any] = {}
    
    @abstractmethod
    def process(self, input_data: Any, **kwargs) -> ProcessingResult:
        """Process input data and return result."""
        pass
    
    @abstractmethod
    def validate_input(self, input_data: Any) -> bool:
        """Validate input data."""
        pass
    
    def run(self, input_data: Any, **kwargs) -> ProcessingResult:
        """
        Main entry point for processing with comprehensive error handling.
        
        Args:
            input_data: Data to process
            **kwargs: Additional processing parameters
            
        Returns:
            ProcessingResult object
        """
        self._start_timer()
        result = ProcessingResult(success=True)
        
        try:
            # Validate input
            if not self.validate_input(input_data):
                result.add_error("Input validation failed")
                return result
            
            # Process data
            self.logger.info(f"Starting processing", input_type=type(input_data).__name__)
            
            processing_result = self.process(input_data, **kwargs)
            result = result.merge(processing_result)
            
            if result.success:
                self.logger.info("Processing completed successfully")
                self._stats["files_processed"] += 1
            else:
                self.logger.error("Processing failed", errors=result.errors)
                self._stats["errors_encountered"] += 1
        
        except Exception as e:
            result.add_error(f"Unexpected error: {str(e)}")
            self.logger.error(f"Unexpected error in processing", error=str(e))
            self._stats["errors_encountered"] += 1
        
        finally:
            processing_time = self._stop_timer()
            result.processing_time = processing_time
            self._stats["total_processing_time"] += processing_time
        
        return result
    
    def process_batch(self, input_items: List[Any], **kwargs) -> List[ProcessingResult]:
        """
        Process multiple items in parallel.
        
        Args:
            input_items: List of items to process
            **kwargs: Additional processing parameters
            
        Returns:
            List of ProcessingResult objects
        """
        results = []
        
        if self.config.max_workers == 1:
            # Sequential processing
            for item in input_items:
                result = self.run(item, **kwargs)
                results.append(result)
        else:
            # Parallel processing
            with ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
                future_to_item = {
                    executor.submit(self.run, item, **kwargs): item 
                    for item in input_items
                }
                
                for future in as_completed(future_to_item):
                    item = future_to_item[future]
                    try:
                        result = future.result(timeout=self.config.timeout)
                        results.append(result)
                    except Exception as e:
                        error_result = ProcessingResult(success=False)
                        error_result.add_error(f"Processing failed for {item}: {str(e)}")
                        results.append(error_result)
        
        return results
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get processing statistics."""
        return {
            **self._stats,
            "processor_name": self.name,
            "cache_size": len(self._cache),
        }
    
    def clear_cache(self) -> None:
        """Clear processing cache."""
        self._cache.clear()
        self.logger.info("Cache cleared")
    
    def save_results(self, data: Any, filename: str, format: str = "json") -> Path:
        """
        Save results to file.
        
        Args:
            data: Data to save
            filename: Output filename
            format: Output format (json, csv, pickle)
            
        Returns:
            Path to saved file
        """
        output_path = self.config.output_dir / filename
        
        try:
            if format == "json":
                with open(output_path, 'w') as f:
                    json.dump(data, f, indent=2, default=str)
            elif format == "csv" and isinstance(data, (list, dict)):
                df = pd.DataFrame(data)
                df.to_csv(output_path, index=False)
            elif format == "pickle":
                with open(output_path, 'wb') as f:
                    pickle.dump(data, f)
            else:
                raise ValueError(f"Unsupported format: {format}")
            
            self.logger.info(f"Results saved", path=str(output_path), format=format)
            return output_path
            
        except Exception as e:
            self.logger.error(f"Failed to save results", error=str(e))
            raise
    
    def _start_timer(self) -> None:
        """Start performance timer."""
        self._start_time = time.time()
    
    def _stop_timer(self) -> float:
        """Stop performance timer and return elapsed time."""
        if self._start_time is None:
            return 0.0
        
        elapsed = time.time() - self._start_time
        self._start_time = None
        return elapsed
    
    @contextmanager
    def temporary_directory(self):
        """Context manager for temporary directory operations."""
        temp_dir = self.config.temp_dir / f"{self.name}_{int(time.time())}"
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            yield temp_dir
        finally:
            if temp_dir.exists():
                shutil.rmtree(temp_dir)


class SequenceProcessor(BaseProcessor):
    """Base class for sequence processing operations."""
    
    def __init__(self, config: BioinformaticsConfig, **kwargs):
        """Initialize sequence processor."""
        super().__init__(config, **kwargs)
        self.sequence_cache: Dict[str, SeqRecord] = {}
    
    def load_fasta(self, fasta_path: Path) -> Dict[str, SeqRecord]:
        """
        Load FASTA file with caching.
        
        Args:
            fasta_path: Path to FASTA file
            
        Returns:
            Dictionary mapping sequence IDs to SeqRecord objects
        """
        cache_key = str(fasta_path)
        
        if self.config.cache_enabled and cache_key in self._cache:
            self.logger.debug(f"Loading FASTA from cache", path=str(fasta_path))
            return self._cache[cache_key]
        
        try:
            sequences = {}
            for record in SeqIO.parse(fasta_path, "fasta"):
                sequences[record.id] = record
            
            if self.config.cache_enabled:
                self._cache[cache_key] = sequences
            
            self.logger.info(f"Loaded FASTA file", 
                           path=str(fasta_path), 
                           sequence_count=len(sequences))
            
            return sequences
            
        except Exception as e:
            self.logger.error(f"Failed to load FASTA file", 
                            path=str(fasta_path), 
                            error=str(e))
            raise
    
    def save_fasta(self, sequences: Dict[str, SeqRecord], output_path: Path) -> None:
        """
        Save sequences to FASTA file.
        
        Args:
            sequences: Dictionary of sequences
            output_path: Output file path
        """
        try:
            with open(output_path, 'w') as f:
                SeqIO.write(sequences.values(), f, "fasta")
            
            self.logger.info(f"Saved FASTA file", 
                           path=str(output_path), 
                           sequence_count=len(sequences))
            
        except Exception as e:
            self.logger.error(f"Failed to save FASTA file", 
                            path=str(output_path), 
                            error=str(e))
            raise
    
    def validate_input(self, input_data: Any) -> bool:
        """Validate sequence input data."""
        if isinstance(input_data, (str, Path)):
            path = Path(input_data)
            return path.exists() and path.suffix.lower() in ['.fa', '.fasta', '.fas']
        elif isinstance(input_data, dict):
            return all(isinstance(v, SeqRecord) for v in input_data.values())
        elif isinstance(input_data, list):
            return all(isinstance(item, SeqRecord) for item in input_data)
        
        return False


class DatabaseProcessor(BaseProcessor):
    """Base class for database operations."""
    
    def __init__(self, config: BioinformaticsConfig, **kwargs):
        """Initialize database processor."""
        super().__init__(config, **kwargs)
        self.connection_pool: Dict[str, Any] = {}
    
    def load_csv_data(self, csv_path: Path, **kwargs) -> pd.DataFrame:
        """
        Load CSV data with caching.
        
        Args:
            csv_path: Path to CSV file
            **kwargs: Additional pandas.read_csv parameters
            
        Returns:
            DataFrame
        """
        cache_key = f"csv_{csv_path}"
        
        if self.config.cache_enabled and cache_key in self._cache:
            self.logger.debug(f"Loading CSV from cache", path=str(csv_path))
            return self._cache[cache_key]
        
        try:
            df = pd.read_csv(csv_path, **kwargs)
            
            if self.config.cache_enabled:
                self._cache[cache_key] = df
            
            self.logger.info(f"Loaded CSV file", 
                           path=str(csv_path), 
                           rows=len(df), 
                           columns=len(df.columns))
            
            return df
            
        except Exception as e:
            self.logger.error(f"Failed to load CSV file", 
                            path=str(csv_path), 
                            error=str(e))
            raise
    
    def save_csv_data(self, df: pd.DataFrame, output_path: Path, **kwargs) -> None:
        """
        Save DataFrame to CSV.
        
        Args:
            df: DataFrame to save
            output_path: Output file path
            **kwargs: Additional pandas.to_csv parameters
        """
        try:
            df.to_csv(output_path, index=False, **kwargs)
            
            self.logger.info(f"Saved CSV file", 
                           path=str(output_path), 
                           rows=len(df), 
                           columns=len(df.columns))
            
        except Exception as e:
            self.logger.error(f"Failed to save CSV file", 
                            path=str(output_path), 
                            error=str(e))
            raise
    
    def validate_input(self, input_data: Any) -> bool:
        """Validate database input data."""
        if isinstance(input_data, (str, Path)):
            path = Path(input_data)
            return path.exists() and path.suffix.lower() in ['.csv', '.tsv']
        elif isinstance(input_data, pd.DataFrame):
            return not input_data.empty
        
        return False


class ParallelProcessor(BaseProcessor):
    """Base class for parallel processing operations."""
    
    def __init__(self, config: BioinformaticsConfig, use_processes: bool = False, **kwargs):
        """
        Initialize parallel processor.
        
        Args:
            config: Processing configuration
            use_processes: Whether to use processes instead of threads
            **kwargs: Additional parameters
        """
        super().__init__(config, **kwargs)
        self.use_processes = use_processes
    
    def execute_parallel(self, 
                        func: callable, 
                        items: List[Any], 
                        **kwargs) -> List[ProcessingResult]:
        """
        Execute function in parallel on items.
        
        Args:
            func: Function to execute
            items: Items to process
            **kwargs: Additional function parameters
            
        Returns:
            List of results
        """
        if self.config.max_workers == 1:
            return [self._execute_with_error_handling(func, item, **kwargs) 
                   for item in items]
        
        executor_class = ProcessPoolExecutor if self.use_processes else ThreadPoolExecutor
        
        with executor_class(max_workers=self.config.max_workers) as executor:
            future_to_item = {
                executor.submit(self._execute_with_error_handling, func, item, **kwargs): item
                for item in items
            }
            
            results = []
            for future in as_completed(future_to_item):
                try:
                    result = future.result(timeout=self.config.timeout)
                    results.append(result)
                except Exception as e:
                    error_result = ProcessingResult(success=False)
                    error_result.add_error(f"Execution failed: {str(e)}")
                    results.append(error_result)
        
        return results
    
    def _execute_with_error_handling(self, 
                                   func: callable, 
                                   item: Any, 
                                   **kwargs) -> ProcessingResult:
        """Execute function with comprehensive error handling."""
        result = ProcessingResult(success=True)
        
        try:
            start_time = time.time()
            data = func(item, **kwargs)
            processing_time = time.time() - start_time
            
            result.data = data
            result.processing_time = processing_time
            
        except Exception as e:
            result.add_error(f"Function execution failed: {str(e)}")
        
        return result
    
    def validate_input(self, input_data: Any) -> bool:
        """Validate parallel processing input."""
        return isinstance(input_data, (list, tuple))


class FileSystemProcessor(BaseProcessor):
    """Base class for file system operations."""
    
    def __init__(self, config: BioinformaticsConfig, **kwargs):
        """Initialize file system processor."""
        super().__init__(config, **kwargs)
    
    def find_files(self, 
                  directory: Path, 
                  pattern: str = "*", 
                  recursive: bool = True) -> List[Path]:
        """
        Find files matching pattern.
        
        Args:
            directory: Directory to search
            pattern: File pattern
            recursive: Whether to search recursively
            
        Returns:
            List of matching file paths
        """
        try:
            if recursive:
                files = list(directory.rglob(pattern))
            else:
                files = list(directory.glob(pattern))
            
            files = [f for f in files if f.is_file()]
            
            self.logger.info(f"Found files", 
                           directory=str(directory), 
                           pattern=pattern, 
                           count=len(files))
            
            return files
            
        except Exception as e:
            self.logger.error(f"Failed to find files", 
                            directory=str(directory), 
                            error=str(e))
            raise
    
    def ensure_directory(self, directory: Path) -> None:
        """Ensure directory exists."""
        try:
            directory.mkdir(parents=True, exist_ok=True)
            self.logger.debug(f"Directory ensured", path=str(directory))
        except Exception as e:
            self.logger.error(f"Failed to create directory", 
                            path=str(directory), 
                            error=str(e))
            raise
    
    def validate_input(self, input_data: Any) -> bool:
        """Validate file system input."""
        if isinstance(input_data, (str, Path)):
            return Path(input_data).exists()
        elif isinstance(input_data, list):
            return all(Path(item).exists() for item in input_data)
        
        return False