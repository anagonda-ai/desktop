"""
Custom exceptions for Plant MGC Analysis Pipeline.

This module defines the exception hierarchy used throughout the application
for better error handling and debugging.
"""

from typing import Optional, Any, Dict


class PlantMGCError(Exception):
    """Base exception class for all Plant MGC Analysis errors."""
    
    def __init__(
        self,
        message: str,
        error_code: Optional[str] = None,
        details: Optional[Dict[str, Any]] = None,
    ):
        """
        Initialize PlantMGCError.
        
        Args:
            message: Human-readable error message
            error_code: Optional error code for programmatic handling
            details: Optional additional error details
        """
        super().__init__(message)
        self.message = message
        self.error_code = error_code
        self.details = details or {}
    
    def __str__(self) -> str:
        """Return string representation of the error."""
        result = self.message
        if self.error_code:
            result = f"[{self.error_code}] {result}"
        if self.details:
            result += f" (Details: {self.details})"
        return result


class AnalysisError(PlantMGCError):
    """Raised when an analysis operation fails."""
    
    def __init__(
        self,
        message: str,
        analysis_type: Optional[str] = None,
        input_data: Optional[str] = None,
        **kwargs
    ):
        """
        Initialize AnalysisError.
        
        Args:
            message: Error message
            analysis_type: Type of analysis that failed
            input_data: Description of input data
            **kwargs: Additional arguments passed to parent
        """
        details = kwargs.get("details", {})
        if analysis_type:
            details["analysis_type"] = analysis_type
        if input_data:
            details["input_data"] = input_data
        
        super().__init__(message, **kwargs)
        self.analysis_type = analysis_type
        self.input_data = input_data


class ConfigurationError(PlantMGCError):
    """Raised when there's a configuration-related error."""
    
    def __init__(
        self,
        message: str,
        config_key: Optional[str] = None,
        config_value: Optional[Any] = None,
        **kwargs
    ):
        """
        Initialize ConfigurationError.
        
        Args:
            message: Error message
            config_key: Configuration key that caused the error
            config_value: Invalid configuration value
            **kwargs: Additional arguments passed to parent
        """
        details = kwargs.get("details", {})
        if config_key:
            details["config_key"] = config_key
        if config_value is not None:
            details["config_value"] = str(config_value)
        
        super().__init__(message, **kwargs)
        self.config_key = config_key
        self.config_value = config_value


class DataError(PlantMGCError):
    """Raised when there's a data-related error."""
    
    def __init__(
        self,
        message: str,
        data_type: Optional[str] = None,
        data_source: Optional[str] = None,
        **kwargs
    ):
        """
        Initialize DataError.
        
        Args:
            message: Error message
            data_type: Type of data that caused the error
            data_source: Source of the problematic data
            **kwargs: Additional arguments passed to parent
        """
        details = kwargs.get("details", {})
        if data_type:
            details["data_type"] = data_type
        if data_source:
            details["data_source"] = data_source
        
        super().__init__(message, **kwargs)
        self.data_type = data_type
        self.data_source = data_source


class ValidationError(PlantMGCError):
    """Raised when input validation fails."""
    
    def __init__(
        self,
        message: str,
        field_name: Optional[str] = None,
        field_value: Optional[Any] = None,
        expected_type: Optional[type] = None,
        **kwargs
    ):
        """
        Initialize ValidationError.
        
        Args:
            message: Error message
            field_name: Name of the field that failed validation
            field_value: Value that failed validation
            expected_type: Expected type for the field
            **kwargs: Additional arguments passed to parent
        """
        details = kwargs.get("details", {})
        if field_name:
            details["field_name"] = field_name
        if field_value is not None:
            details["field_value"] = str(field_value)
        if expected_type:
            details["expected_type"] = expected_type.__name__
        
        super().__init__(message, **kwargs)
        self.field_name = field_name
        self.field_value = field_value
        self.expected_type = expected_type


class DatabaseError(PlantMGCError):
    """Raised when database operations fail."""
    
    def __init__(
        self,
        message: str,
        database_name: Optional[str] = None,
        operation: Optional[str] = None,
        **kwargs
    ):
        """
        Initialize DatabaseError.
        
        Args:
            message: Error message
            database_name: Name of the database
            operation: Database operation that failed
            **kwargs: Additional arguments passed to parent
        """
        details = kwargs.get("details", {})
        if database_name:
            details["database_name"] = database_name
        if operation:
            details["operation"] = operation
        
        super().__init__(message, **kwargs)
        self.database_name = database_name
        self.operation = operation


class NetworkError(PlantMGCError):
    """Raised when network operations fail."""
    
    def __init__(
        self,
        message: str,
        url: Optional[str] = None,
        status_code: Optional[int] = None,
        **kwargs
    ):
        """
        Initialize NetworkError.
        
        Args:
            message: Error message
            url: URL that caused the error
            status_code: HTTP status code
            **kwargs: Additional arguments passed to parent
        """
        details = kwargs.get("details", {})
        if url:
            details["url"] = url
        if status_code:
            details["status_code"] = status_code
        
        super().__init__(message, **kwargs)
        self.url = url
        self.status_code = status_code


class ComputeError(PlantMGCError):
    """Raised when computational operations fail."""
    
    def __init__(
        self,
        message: str,
        computation_type: Optional[str] = None,
        resource_type: Optional[str] = None,
        **kwargs
    ):
        """
        Initialize ComputeError.
        
        Args:
            message: Error message
            computation_type: Type of computation that failed
            resource_type: Type of resource that was unavailable
            **kwargs: Additional arguments passed to parent
        """
        details = kwargs.get("details", {})
        if computation_type:
            details["computation_type"] = computation_type
        if resource_type:
            details["resource_type"] = resource_type
        
        super().__init__(message, **kwargs)
        self.computation_type = computation_type
        self.resource_type = resource_type


class FileSystemError(PlantMGCError):
    """Raised when file system operations fail."""
    
    def __init__(
        self,
        message: str,
        file_path: Optional[str] = None,
        operation: Optional[str] = None,
        **kwargs
    ):
        """
        Initialize FileSystemError.
        
        Args:
            message: Error message
            file_path: Path to the file that caused the error
            operation: File operation that failed
            **kwargs: Additional arguments passed to parent
        """
        details = kwargs.get("details", {})
        if file_path:
            details["file_path"] = file_path
        if operation:
            details["operation"] = operation
        
        super().__init__(message, **kwargs)
        self.file_path = file_path
        self.operation = operation


class ProcessingError(PlantMGCError):
    """Raised when general processing operations fail."""
    
    def __init__(
        self,
        message: str,
        processing_type: Optional[str] = None,
        **kwargs
    ):
        """
        Initialize ProcessingError.
        
        Args:
            message: Error message
            processing_type: Type of processing that failed
            **kwargs: Additional arguments passed to parent
        """
        details = kwargs.get("details", {})
        if processing_type:
            details["processing_type"] = processing_type
        
        super().__init__(message, **kwargs)
        self.processing_type = processing_type