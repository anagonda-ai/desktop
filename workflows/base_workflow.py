"""
Base workflow class providing unified interface for all analysis workflows.

All bioinformatics workflows inherit from BaseWorkflow to ensure consistent
behavior, error handling, and result management across the system.
"""

import time
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Any, List, Optional, Union
from datetime import datetime

from core import Config, get_config
from core.exceptions import BioinformaticsError, ProcessingError
from core.types import ProcessingStatus
from utils import FileManager, ProgressTracker


@dataclass
class WorkflowResult:
    """Standardized result object for all workflows."""
    workflow_name: str
    status: ProcessingStatus
    start_time: datetime
    end_time: Optional[datetime] = None
    execution_time: Optional[float] = None
    results: Dict[str, Any] = field(default_factory=dict)
    metrics: Dict[str, Any] = field(default_factory=dict)
    artifacts: List[Path] = field(default_factory=list)
    error_message: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert result to dictionary."""
        return {
            'workflow_name': self.workflow_name,
            'status': self.status.value,
            'start_time': self.start_time.isoformat(),
            'end_time': self.end_time.isoformat() if self.end_time else None,
            'execution_time': self.execution_time,
            'results': self.results,
            'metrics': self.metrics,
            'artifacts': [str(p) for p in self.artifacts],
            'error_message': self.error_message
        }


class BaseWorkflow(ABC):
    """
    Base class for all bioinformatics analysis workflows.
    
    Provides unified interface, error handling, progress tracking,
    and result management for all analysis pipelines.
    """
    
    def __init__(self, name: str):
        self.name = name
        self.config = get_config()
        self.file_manager = FileManager()
        self.logger = logging.getLogger(f"workflows.{name}")
        
        # Execution state
        self.current_result: Optional[WorkflowResult] = None
        self.workspace: Optional[Path] = None
        
    def setup_workspace(self, workspace: Path) -> None:
        """Setup workflow-specific workspace."""
        self.workspace = workspace
        workflow_dir = workspace / "results" / self.name.lower()
        self.file_manager.ensure_directory(workflow_dir)
        
        # Create workflow-specific subdirectories
        for subdir in self.get_required_subdirectories():
            self.file_manager.ensure_directory(workflow_dir / subdir)
    
    def get_required_subdirectories(self) -> List[str]:
        """Get list of required subdirectories for this workflow."""
        return ["raw", "processed", "final", "temp"]
    
    def validate_parameters(self, params: Dict[str, Any]) -> None:
        """Validate workflow parameters."""
        required_params = self.get_required_parameters()
        missing_params = [p for p in required_params if p not in params]
        
        if missing_params:
            raise BioinformaticsError(
                f"Missing required parameters for {self.name} workflow: {missing_params}"
            )
    
    @abstractmethod
    def get_required_parameters(self) -> List[str]:
        """Get list of required parameters for this workflow."""
        pass
    
    @abstractmethod
    def _execute_workflow(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Execute the workflow logic. Must be implemented by subclasses."""
        pass
    
    def pre_execution_hook(self, params: Dict[str, Any]) -> None:
        """Hook called before workflow execution."""
        pass
    
    def post_execution_hook(self, result: WorkflowResult) -> None:
        """Hook called after workflow execution."""
        pass
    
    def cleanup_temporary_files(self) -> None:
        """Clean up temporary files created during execution."""
        if self.workspace:
            temp_dir = self.workspace / "results" / self.name.lower() / "temp"
            if temp_dir.exists():
                import shutil
                shutil.rmtree(temp_dir)
                self.logger.debug(f"Cleaned up temporary files in {temp_dir}")
    
    def save_workflow_metadata(self, result: WorkflowResult) -> Path:
        """Save workflow execution metadata."""
        if not self.workspace:
            raise ProcessingError("Workspace not initialized")
        
        metadata_file = (
            self.workspace / "results" / self.name.lower() / 
            f"metadata_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        )
        
        import json
        with open(metadata_file, 'w') as f:
            json.dump(result.to_dict(), f, indent=2)
        
        self.logger.info(f"Saved workflow metadata to {metadata_file}")
        return metadata_file
    
    def get_output_paths(self) -> Dict[str, Path]:
        """Get standardized output paths for this workflow."""
        if not self.workspace:
            raise ProcessingError("Workspace not initialized")
        
        base_dir = self.workspace / "results" / self.name.lower()
        
        return {
            'base': base_dir,
            'raw': base_dir / "raw",
            'processed': base_dir / "processed", 
            'final': base_dir / "final",
            'temp': base_dir / "temp"
        }
    
    def collect_artifacts(self) -> List[Path]:
        """Collect all output artifacts produced by this workflow."""
        artifacts = []
        output_paths = self.get_output_paths()
        
        for path_type, path in output_paths.items():
            if path_type != 'temp' and path.exists():
                # Add all files in the directory
                for file_path in path.rglob('*'):
                    if file_path.is_file():
                        artifacts.append(file_path)
        
        return artifacts
    
    def calculate_metrics(self, execution_time: float) -> Dict[str, Any]:
        """Calculate workflow execution metrics."""
        artifacts = self.collect_artifacts()
        
        # Calculate total output size
        total_size = sum(f.stat().st_size for f in artifacts if f.exists())
        
        return {
            'execution_time_seconds': execution_time,
            'total_artifacts': len(artifacts),
            'total_output_size_bytes': total_size,
            'artifacts_by_type': self._categorize_artifacts(artifacts)
        }
    
    def _categorize_artifacts(self, artifacts: List[Path]) -> Dict[str, int]:
        """Categorize artifacts by file type."""
        categories = {}
        for artifact in artifacts:
            ext = artifact.suffix.lower()
            categories[ext] = categories.get(ext, 0) + 1
        return categories
    
    def execute(self, params: Dict[str, Any]) -> WorkflowResult:
        """
        Execute the workflow with unified error handling and result management.
        
        Args:
            params: Workflow parameters
            
        Returns:
            WorkflowResult object with execution details
        """
        start_time = datetime.now()
        self.logger.info(f"Starting {self.name} workflow")
        
        # Initialize result
        self.current_result = WorkflowResult(
            workflow_name=self.name,
            status=ProcessingStatus.RUNNING,
            start_time=start_time
        )
        
        try:
            # Setup workspace if provided
            if 'workspace' in params:
                self.setup_workspace(Path(params['workspace']))
            
            # Validate parameters
            self.validate_parameters(params)
            
            # Pre-execution hook
            self.pre_execution_hook(params)
            
            # Execute workflow
            start_exec_time = time.time()
            results = self._execute_workflow(params)
            execution_time = time.time() - start_exec_time
            
            # Update result
            end_time = datetime.now()
            self.current_result.status = ProcessingStatus.COMPLETED
            self.current_result.end_time = end_time
            self.current_result.execution_time = execution_time
            self.current_result.results = results
            self.current_result.metrics = self.calculate_metrics(execution_time)
            self.current_result.artifacts = self.collect_artifacts()
            
            # Post-execution hook
            self.post_execution_hook(self.current_result)
            
            # Save metadata
            if self.workspace:
                metadata_file = self.save_workflow_metadata(self.current_result)
                self.current_result.artifacts.append(metadata_file)
            
            # Cleanup
            self.cleanup_temporary_files()
            
            self.logger.info(
                f"Completed {self.name} workflow in {execution_time:.2f}s "
                f"with {len(self.current_result.artifacts)} artifacts"
            )
            
            return self.current_result
            
        except Exception as e:
            # Handle execution failure
            end_time = datetime.now()
            execution_time = (end_time - start_time).total_seconds()
            
            self.current_result.status = ProcessingStatus.FAILED
            self.current_result.end_time = end_time
            self.current_result.execution_time = execution_time
            self.current_result.error_message = str(e)
            
            self.logger.error(f"Workflow {self.name} failed: {e}")
            
            # Try to save error metadata
            if self.workspace:
                try:
                    self.save_workflow_metadata(self.current_result)
                except Exception as meta_error:
                    self.logger.error(f"Failed to save error metadata: {meta_error}")
            
            # Cleanup on failure
            self.cleanup_temporary_files()
            
            raise ProcessingError(
                f"Workflow {self.name} failed: {e}",
                processor_name=self.name,
                stage="workflow_execution",
                original_error=e
            )
    
    def get_status(self) -> Optional[ProcessingStatus]:
        """Get current workflow status."""
        return self.current_result.status if self.current_result else None
    
    def get_progress(self) -> Dict[str, Any]:
        """Get current workflow progress information."""
        if not self.current_result:
            return {"status": "not_started"}
        
        progress = {
            "status": self.current_result.status.value,
            "start_time": self.current_result.start_time.isoformat(),
            "workflow_name": self.name
        }
        
        if self.current_result.execution_time:
            progress["execution_time"] = self.current_result.execution_time
        
        if self.current_result.error_message:
            progress["error_message"] = self.current_result.error_message
        
        return progress