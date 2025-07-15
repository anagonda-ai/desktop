"""
Unified logging configuration for bioinformatics applications.

Provides centralized logging setup with file rotation, structured output,
and performance monitoring capabilities.
"""

import sys
import logging
import logging.handlers
from pathlib import Path
from typing import Optional, Dict, Any
from datetime import datetime
import json


class StructuredFormatter(logging.Formatter):
    """
    Structured logging formatter that outputs JSON format.
    
    Provides structured logs with metadata for better analysis
    and monitoring capabilities.
    """
    
    def format(self, record: logging.LogRecord) -> str:
        """Format log record as structured JSON."""
        log_entry = {
            'timestamp': datetime.fromtimestamp(record.created).isoformat(),
            'level': record.levelname,
            'logger': record.name,
            'message': record.getMessage(),
            'module': record.module,
            'function': record.funcName,
            'line': record.lineno
        }
        
        # Add exception information if present
        if record.exc_info:
            log_entry['exception'] = self.formatException(record.exc_info)
        
        # Add extra fields from record
        for key, value in record.__dict__.items():
            if key not in ['name', 'msg', 'args', 'levelname', 'levelno', 'pathname',
                          'filename', 'module', 'lineno', 'funcName', 'created',
                          'msecs', 'relativeCreated', 'thread', 'threadName',
                          'processName', 'process', 'getMessage', 'exc_info',
                          'exc_text', 'stack_info']:
                log_entry[key] = value
        
        return json.dumps(log_entry)


class BioinformaticsLogger:
    """
    Centralized logger for bioinformatics applications.
    
    Provides unified logging configuration with file rotation,
    structured output, and performance monitoring.
    """
    
    def __init__(self, name: str = "bioinformatics"):
        self.name = name
        self.logger = logging.getLogger(name)
        self.handlers_configured = False
    
    def setup_logging(
        self,
        level: str = "INFO",
        log_dir: Optional[Path] = None,
        enable_console: bool = True,
        enable_file: bool = True,
        enable_structured: bool = False,
        max_file_size: int = 10 * 1024 * 1024,  # 10MB
        backup_count: int = 5
    ) -> None:
        """
        Setup logging configuration.
        
        Args:
            level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
            log_dir: Directory for log files
            enable_console: Enable console output
            enable_file: Enable file output
            enable_structured: Enable structured JSON logging
            max_file_size: Maximum file size before rotation
            backup_count: Number of backup files to keep
        """
        if self.handlers_configured:
            return
        
        # Clear any existing handlers
        self.logger.handlers.clear()
        
        # Set logging level
        log_level = getattr(logging, level.upper())
        self.logger.setLevel(log_level)
        
        # Create formatters
        console_formatter = logging.Formatter(
            fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        file_formatter = logging.Formatter(
            fmt='%(asctime)s - %(name)s - %(levelname)s - %(module)s:%(funcName)s:%(lineno)d - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        structured_formatter = StructuredFormatter()
        
        # Console handler
        if enable_console:
            console_handler = logging.StreamHandler(sys.stdout)
            console_handler.setLevel(log_level)
            
            if enable_structured:
                console_handler.setFormatter(structured_formatter)
            else:
                console_handler.setFormatter(console_formatter)
            
            self.logger.addHandler(console_handler)
        
        # File handlers
        if enable_file and log_dir:
            log_dir = Path(log_dir)
            log_dir.mkdir(parents=True, exist_ok=True)
            
            # Regular log file
            log_file = log_dir / f"{self.name}.log"
            file_handler = logging.handlers.RotatingFileHandler(
                log_file,
                maxBytes=max_file_size,
                backupCount=backup_count
            )
            file_handler.setLevel(log_level)
            file_handler.setFormatter(file_formatter)
            self.logger.addHandler(file_handler)
            
            # Error log file
            error_file = log_dir / f"{self.name}_errors.log"
            error_handler = logging.handlers.RotatingFileHandler(
                error_file,
                maxBytes=max_file_size,
                backupCount=backup_count
            )
            error_handler.setLevel(logging.ERROR)
            error_handler.setFormatter(file_formatter)
            self.logger.addHandler(error_handler)
            
            # Structured log file (if enabled)
            if enable_structured:
                structured_file = log_dir / f"{self.name}_structured.jsonl"
                structured_handler = logging.handlers.RotatingFileHandler(
                    structured_file,
                    maxBytes=max_file_size,
                    backupCount=backup_count
                )
                structured_handler.setLevel(log_level)
                structured_handler.setFormatter(structured_formatter)
                self.logger.addHandler(structured_handler)
        
        self.handlers_configured = True
        
        # Log setup completion
        self.logger.info(f"Logging configured - Level: {level}, Console: {enable_console}, File: {enable_file}")
    
    def get_logger(self, name: Optional[str] = None) -> logging.Logger:
        """Get a logger instance."""
        if name:
            return logging.getLogger(f"{self.name}.{name}")
        return self.logger
    
    def log_performance(self, operation: str, duration: float, metadata: Optional[Dict[str, Any]] = None):
        """Log performance metrics."""
        perf_logger = self.get_logger("performance")
        
        message = f"PERFORMANCE: {operation} completed in {duration:.3f}s"
        
        if metadata:
            # Add metadata as extra fields
            perf_logger.info(message, extra=metadata)
        else:
            perf_logger.info(message)
    
    def log_workflow_start(self, workflow_name: str, params: Dict[str, Any]):
        """Log workflow start."""
        workflow_logger = self.get_logger("workflow")
        workflow_logger.info(
            f"Starting {workflow_name} workflow",
            extra={'workflow': workflow_name, 'parameters': params}
        )
    
    def log_workflow_end(self, workflow_name: str, duration: float, status: str, results: Optional[Dict[str, Any]] = None):
        """Log workflow completion."""
        workflow_logger = self.get_logger("workflow")
        
        extra_data = {
            'workflow': workflow_name,
            'duration': duration,
            'status': status
        }
        
        if results:
            extra_data['results_summary'] = results
        
        workflow_logger.info(
            f"Completed {workflow_name} workflow - Status: {status}, Duration: {duration:.2f}s",
            extra=extra_data
        )
    
    def log_error_with_context(self, error: Exception, context: Dict[str, Any]):
        """Log error with contextual information."""
        error_logger = self.get_logger("error")
        error_logger.error(
            f"Error occurred: {error}",
            extra={'error_type': type(error).__name__, 'context': context},
            exc_info=True
        )


# Global logger instance
_global_logger: Optional[BioinformaticsLogger] = None


def setup_logging(
    level: str = "INFO",
    log_dir: Optional[Path] = None,
    enable_console: bool = True,
    enable_file: bool = True,
    enable_structured: bool = False,
    logger_name: str = "bioinformatics"
) -> BioinformaticsLogger:
    """
    Setup global logging configuration.
    
    Args:
        level: Logging level
        log_dir: Directory for log files
        enable_console: Enable console output
        enable_file: Enable file output
        enable_structured: Enable structured JSON logging
        logger_name: Name for the logger
        
    Returns:
        Configured BioinformaticsLogger instance
    """
    global _global_logger
    
    if _global_logger is None:
        _global_logger = BioinformaticsLogger(logger_name)
    
    _global_logger.setup_logging(
        level=level,
        log_dir=log_dir,
        enable_console=enable_console,
        enable_file=enable_file,
        enable_structured=enable_structured
    )
    
    return _global_logger


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """
    Get a logger instance.
    
    Args:
        name: Optional logger name suffix
        
    Returns:
        Logger instance
    """
    global _global_logger
    
    if _global_logger is None:
        _global_logger = setup_logging()
    
    return _global_logger.get_logger(name)


def log_performance(operation: str, duration: float, metadata: Optional[Dict[str, Any]] = None):
    """Log performance metrics."""
    global _global_logger
    
    if _global_logger is None:
        _global_logger = setup_logging()
    
    _global_logger.log_performance(operation, duration, metadata)


def log_workflow_start(workflow_name: str, params: Dict[str, Any]):
    """Log workflow start."""
    global _global_logger
    
    if _global_logger is None:
        _global_logger = setup_logging()
    
    _global_logger.log_workflow_start(workflow_name, params)


def log_workflow_end(workflow_name: str, duration: float, status: str, results: Optional[Dict[str, Any]] = None):
    """Log workflow completion."""
    global _global_logger
    
    if _global_logger is None:
        _global_logger = setup_logging()
    
    _global_logger.log_workflow_end(workflow_name, duration, status, results)