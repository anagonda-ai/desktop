"""
Logging utilities for Plant MGC Analysis Pipeline.

This module provides centralized logging configuration using loguru
with structured logging support and performance monitoring.
"""

import sys
from pathlib import Path
from typing import Optional, Dict, Any
from loguru import logger


def setup_logging(
    level: str = "INFO",
    log_file: Optional[Path] = None,
    format_string: Optional[str] = None,
    enable_json: bool = False,
    rotation: str = "1 week",
    retention: str = "1 month",
    compression: str = "gz",
) -> None:
    """
    Setup centralized logging configuration.
    
    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Optional path to log file
        format_string: Custom format string
        enable_json: Enable JSON structured logging
        rotation: Log rotation interval
        retention: Log retention period
        compression: Compression format for rotated logs
    """
    # Remove default logger
    logger.remove()
    
    # Default format
    if format_string is None:
        if enable_json:
            format_string = "{message}"
        else:
            format_string = (
                "{time:YYYY-MM-DD HH:mm:ss.SSS} | "
                "{level: <8} | "
                "{name}:{function}:{line} | "
                "{message}"
            )
    
    # Console handler
    logger.add(
        sys.stderr,
        format=format_string,
        level=level,
        colorize=True,
        serialize=enable_json,
    )
    
    # File handler
    if log_file:
        logger.add(
            log_file,
            format=format_string,
            level=level,
            rotation=rotation,
            retention=retention,
            compression=compression,
            serialize=enable_json,
        )
    
    # Set up error file handler
    if log_file:
        error_file = log_file.parent / f"{log_file.stem}_errors.log"
        logger.add(
            error_file,
            format=format_string,
            level="ERROR",
            rotation=rotation,
            retention=retention,
            compression=compression,
        )
    
    logger.info("Logging system initialized")


def log_function_call(func_name: str, **kwargs: Any) -> None:
    """
    Log function call with parameters.
    
    Args:
        func_name: Name of the function
        **kwargs: Function parameters
    """
    logger.debug(f"Calling {func_name} with parameters: {kwargs}")


def log_performance(func_name: str, execution_time: float, **metrics: Any) -> None:
    """
    Log performance metrics.
    
    Args:
        func_name: Name of the function
        execution_time: Execution time in seconds
        **metrics: Additional performance metrics
    """
    logger.info(
        f"Performance: {func_name} executed in {execution_time:.3f}s",
        extra={"execution_time": execution_time, "metrics": metrics}
    )


def log_analysis_start(analysis_type: str, input_data: Dict[str, Any]) -> None:
    """
    Log analysis start.
    
    Args:
        analysis_type: Type of analysis
        input_data: Input data description
    """
    logger.info(
        f"Starting {analysis_type} analysis",
        extra={"analysis_type": analysis_type, "input_data": input_data}
    )


def log_analysis_complete(
    analysis_type: str, 
    results: Dict[str, Any], 
    execution_time: float
) -> None:
    """
    Log analysis completion.
    
    Args:
        analysis_type: Type of analysis
        results: Analysis results summary
        execution_time: Execution time in seconds
    """
    logger.info(
        f"Completed {analysis_type} analysis in {execution_time:.3f}s",
        extra={
            "analysis_type": analysis_type,
            "results": results,
            "execution_time": execution_time
        }
    )


def log_error_with_context(
    error: Exception, 
    context: Dict[str, Any], 
    analysis_type: Optional[str] = None
) -> None:
    """
    Log error with additional context.
    
    Args:
        error: Exception that occurred
        context: Additional context information
        analysis_type: Optional analysis type
    """
    logger.error(
        f"Error in {analysis_type or 'operation'}: {error}",
        extra={
            "error_type": type(error).__name__,
            "error_message": str(error),
            "context": context,
            "analysis_type": analysis_type
        }
    )


class LoggerMixin:
    """Mixin class to add logging capabilities to any class."""
    
    @property
    def logger(self):
        """Get logger instance with class name."""
        return logger.bind(class_name=self.__class__.__name__)
    
    def log_method_call(self, method_name: str, **kwargs: Any) -> None:
        """Log method call with parameters."""
        self.logger.debug(f"Calling {method_name} with parameters: {kwargs}")
    
    def log_method_complete(self, method_name: str, execution_time: float) -> None:
        """Log method completion."""
        self.logger.debug(f"Method {method_name} completed in {execution_time:.3f}s")
    
    def log_error(self, error: Exception, method_name: str, **context: Any) -> None:
        """Log error with method context."""
        self.logger.error(
            f"Error in {method_name}: {error}",
            extra={
                "method_name": method_name,
                "error_type": type(error).__name__,
                "context": context
            }
        )


def performance_monitor(func):
    """
    Decorator to monitor function performance.
    
    Args:
        func: Function to monitor
        
    Returns:
        Wrapped function with performance monitoring
    """
    import time
    import functools
    
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        
        try:
            result = func(*args, **kwargs)
            execution_time = time.time() - start_time
            
            log_performance(
                func.__name__,
                execution_time,
                args_count=len(args),
                kwargs_count=len(kwargs)
            )
            
            return result
            
        except Exception as e:
            execution_time = time.time() - start_time
            log_error_with_context(
                e,
                {
                    "function": func.__name__,
                    "args_count": len(args),
                    "kwargs_count": len(kwargs),
                    "execution_time": execution_time
                }
            )
            raise
    
    return wrapper