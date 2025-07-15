"""
Custom exception hierarchy for bioinformatics operations.

This module provides a structured exception system with detailed error context
and error codes for better debugging and error handling.
"""

from typing import Any, Dict, Optional
from enum import Enum


class ErrorCode(Enum):
    """Error codes for different types of failures."""
    # File operation errors
    FILE_NOT_FOUND = "FILE_001"
    FILE_PERMISSION_DENIED = "FILE_002"
    FILE_FORMAT_INVALID = "FILE_003"
    FILE_CORRUPTED = "FILE_004"
    
    # API errors
    API_CONNECTION_FAILED = "API_001"
    API_RATE_LIMIT = "API_002"
    API_AUTHENTICATION = "API_003"
    API_INVALID_RESPONSE = "API_004"
    
    # Validation errors
    INVALID_SEQUENCE = "VAL_001"
    INVALID_COORDINATES = "VAL_002"
    INVALID_PARAMETERS = "VAL_003"
    MISSING_REQUIRED_FIELD = "VAL_004"
    
    # Processing errors
    PROCESSING_FAILED = "PROC_001"
    INSUFFICIENT_MEMORY = "PROC_002"
    TIMEOUT = "PROC_003"
    DEPENDENCY_MISSING = "PROC_004"


class BioinformaticsError(Exception):
    """
    Base exception for all bioinformatics operations.
    
    Provides structured error handling with error codes, context, and suggestions.
    """
    
    def __init__(
        self,
        message: str,
        error_code: Optional[ErrorCode] = None,
        context: Optional[Dict[str, Any]] = None,
        suggestion: Optional[str] = None,
        original_error: Optional[Exception] = None
    ):
        super().__init__(message)
        self.error_code = error_code
        self.context = context or {}
        self.suggestion = suggestion
        self.original_error = original_error
    
    def __str__(self) -> str:
        parts = [super().__str__()]
        
        if self.error_code:
            parts.append(f"Error Code: {self.error_code.value}")
        
        if self.context:
            context_str = ", ".join(f"{k}={v}" for k, v in self.context.items())
            parts.append(f"Context: {context_str}")
        
        if self.suggestion:
            parts.append(f"Suggestion: {self.suggestion}")
        
        return " | ".join(parts)


class ValidationError(BioinformaticsError):
    """Raised when input validation fails."""
    
    def __init__(
        self,
        message: str,
        field_name: Optional[str] = None,
        invalid_value: Any = None,
        **kwargs
    ):
        context = kwargs.get('context', {})
        if field_name:
            context['field'] = field_name
        if invalid_value is not None:
            context['value'] = str(invalid_value)
        
        kwargs['context'] = context
        super().__init__(message, **kwargs)


class FileOperationError(BioinformaticsError):
    """Raised when file operations fail."""
    
    def __init__(
        self,
        message: str,
        file_path: Optional[str] = None,
        operation: Optional[str] = None,
        **kwargs
    ):
        context = kwargs.get('context', {})
        if file_path:
            context['file_path'] = file_path
        if operation:
            context['operation'] = operation
            
        kwargs['context'] = context
        super().__init__(message, **kwargs)


class APIError(BioinformaticsError):
    """Raised when API operations fail."""
    
    def __init__(
        self,
        message: str,
        api_name: Optional[str] = None,
        endpoint: Optional[str] = None,
        status_code: Optional[int] = None,
        **kwargs
    ):
        context = kwargs.get('context', {})
        if api_name:
            context['api'] = api_name
        if endpoint:
            context['endpoint'] = endpoint
        if status_code:
            context['status_code'] = status_code
            
        kwargs['context'] = context
        super().__init__(message, **kwargs)


class ProcessingError(BioinformaticsError):
    """Raised when processing operations fail."""
    
    def __init__(
        self,
        message: str,
        processor_name: Optional[str] = None,
        stage: Optional[str] = None,
        **kwargs
    ):
        context = kwargs.get('context', {})
        if processor_name:
            context['processor'] = processor_name
        if stage:
            context['stage'] = stage
            
        kwargs['context'] = context
        super().__init__(message, **kwargs)


class ConfigurationError(BioinformaticsError):
    """Raised when configuration is invalid or missing."""
    pass


class DependencyError(BioinformaticsError):
    """Raised when required dependencies are missing or invalid."""
    pass