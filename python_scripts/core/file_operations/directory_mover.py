"""
Industry-Level Directory Mover.

This module provides object-oriented directory moving functionality with comprehensive
error handling, validation, backup capabilities, and batch processing.
"""

import os
import shutil
from typing import Dict, List, Optional, Any, Set
from dataclasses import dataclass, field
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from ..base_classes import (
    BaseProcessor,
    BioinformaticsConfig,
    ProcessingResult,
    FileSystemProcessor
)


@dataclass
class DirectoryMoveConfig:
    """Configuration for directory moving operations."""
    
    exclude_patterns: List[str] = field(default_factory=list)
    include_patterns: List[str] = field(default_factory=list)
    create_backup: bool = False
    backup_suffix: str = ".backup"
    dry_run: bool = False
    overwrite_existing: bool = False
    preserve_permissions: bool = True
    verify_moves: bool = True
    max_retries: int = 3


@dataclass
class MoveOperation:
    """Information about a directory move operation."""
    
    source_path: Path
    destination_path: Path
    backup_created: bool = False
    status: str = "pending"  # pending, success, failed, skipped
    error_message: Optional[str] = None
    size_bytes: int = 0
    
    @property
    def source_name(self) -> str:
        """Get source directory name."""
        return self.source_path.name
    
    @property
    def destination_name(self) -> str:
        """Get destination directory name."""
        return self.destination_path.name


class DirectoryMover(FileSystemProcessor):
    """
    Professional directory mover with comprehensive functionality.
    """
    
    def __init__(self, 
                 config: BioinformaticsConfig,
                 move_config: DirectoryMoveConfig):
        """
        Initialize directory mover.
        
        Args:
            config: Base processing configuration
            move_config: Move-specific configuration
        """
        super().__init__(config, name="DirectoryMover")
        
        self.move_config = move_config
        
        # Statistics
        self.move_stats = {
            "directories_found": 0,
            "directories_moved": 0,
            "directories_skipped": 0,
            "directories_failed": 0,
            "total_size_moved_mb": 0.0,
            "backups_created": 0,
        }
    
    def process(self, input_data: Path, **kwargs) -> ProcessingResult:
        """
        Process directory moving operations.
        
        Args:
            input_data: Path to source directory
            **kwargs: Additional processing parameters including:
                - destination_dir: Target directory path
                - exclude_directories: List of directory names to exclude
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            source_dir = Path(input_data) if isinstance(input_data, str) else input_data
            destination_dir = Path(kwargs.get('destination_dir'))
            exclude_dirs = set(kwargs.get('exclude_directories', []))
            
            self.logger.info(f"Starting directory moving", 
                           source=str(source_dir),
                           destination=str(destination_dir),
                           dry_run=self.move_config.dry_run)
            
            # Validate inputs
            if not source_dir.exists():
                result.add_error(f"Source directory does not exist: {source_dir}")
                return result
            
            if not source_dir.is_dir():
                result.add_error(f"Source path is not a directory: {source_dir}")
                return result
            
            # Step 1: Find directories to move
            move_operations = self._plan_move_operations(source_dir, destination_dir, exclude_dirs)
            self.move_stats["directories_found"] = len(move_operations)
            
            if not move_operations:
                result.add_warning("No directories found to move")
                return result
            
            # Step 2: Execute move operations
            execution_result = self._execute_move_operations(move_operations)
            result = result.merge(execution_result)
            
            # Update metadata
            result.metadata.update({
                "move_stats": self.move_stats,
                "source_directory": str(source_dir),
                "destination_directory": str(destination_dir),
                "excluded_directories": list(exclude_dirs),
                "dry_run": self.move_config.dry_run,
                "operations_planned": len(move_operations)
            })
            
            self.logger.info(f"Directory moving completed",
                           directories_moved=self.move_stats["directories_moved"],
                           directories_failed=self.move_stats["directories_failed"])
            
        except Exception as e:
            result.add_error(f"Directory moving failed: {str(e)}")
        
        return result
    
    def _plan_move_operations(self, 
                             source_dir: Path, 
                             destination_dir: Path, 
                             exclude_dirs: Set[str]) -> List[MoveOperation]:
        """Plan move operations for directories."""
        operations = []
        
        try:
            for item_path in source_dir.iterdir():
                if not item_path.is_dir():
                    continue
                
                # Check exclusions
                if item_path.name in exclude_dirs:
                    self.logger.debug(f"Excluding directory", dir=item_path.name)
                    continue
                
                # Check exclude patterns
                if self._should_exclude_directory(item_path):
                    self.logger.debug(f"Excluding directory by pattern", dir=item_path.name)
                    continue
                
                # Check include patterns if specified
                if (self.move_config.include_patterns and 
                    not self._should_include_directory(item_path)):
                    self.logger.debug(f"Not including directory by pattern", dir=item_path.name)
                    continue
                
                # Calculate destination path
                dest_path = destination_dir / item_path.name
                
                # Check for conflicts
                if dest_path.exists() and not self.move_config.overwrite_existing:
                    self.logger.warning(f"Destination already exists", 
                                      source=str(item_path),
                                      destination=str(dest_path))
                    continue
                
                # Calculate directory size
                size_bytes = self._calculate_directory_size(item_path)
                
                operation = MoveOperation(
                    source_path=item_path,
                    destination_path=dest_path,
                    size_bytes=size_bytes
                )
                operations.append(operation)
                
                self.logger.debug(f"Planned move operation", 
                                source=item_path.name,
                                destination=str(dest_path),
                                size_mb=size_bytes / (1024 * 1024))
            
            self.logger.info(f"Planned move operations", count=len(operations))
            return operations
            
        except Exception as e:
            self.logger.error(f"Failed to plan move operations", error=str(e))
            raise
    
    def _should_exclude_directory(self, dir_path: Path) -> bool:
        """Check if directory should be excluded based on patterns."""
        if not self.move_config.exclude_patterns:
            return False
        
        dir_name = dir_path.name
        for pattern in self.move_config.exclude_patterns:
            if pattern in dir_name or dir_name.startswith(pattern) or dir_name.endswith(pattern):
                return True
        
        return False
    
    def _should_include_directory(self, dir_path: Path) -> bool:
        """Check if directory should be included based on patterns."""
        if not self.move_config.include_patterns:
            return True
        
        dir_name = dir_path.name
        for pattern in self.move_config.include_patterns:
            if pattern in dir_name or dir_name.startswith(pattern) or dir_name.endswith(pattern):
                return True
        
        return False
    
    def _calculate_directory_size(self, dir_path: Path) -> int:
        """Calculate total size of directory in bytes."""
        try:
            total_size = 0
            for root, dirs, files in os.walk(dir_path):
                for file in files:
                    file_path = Path(root) / file
                    try:
                        total_size += file_path.stat().st_size
                    except (OSError, FileNotFoundError):
                        # Skip files that can't be accessed
                        pass
            return total_size
        except Exception:
            return 0
    
    def _execute_move_operations(self, operations: List[MoveOperation]) -> ProcessingResult:
        """Execute planned move operations."""
        result = ProcessingResult(success=True)
        
        try:
            if self.move_config.dry_run:
                self.logger.info(f"DRY RUN: Would execute move operations", count=len(operations))
                
                for operation in operations:
                    operation.status = "dry_run"
                    self.move_stats["directories_moved"] += 1
                    
                    self.logger.info(f"DRY RUN: Would move directory", 
                                   source=operation.source_name,
                                   destination=str(operation.destination_path))
            else:
                # Create destination directory
                destination_parent = operations[0].destination_path.parent if operations else None
                if destination_parent:
                    destination_parent.mkdir(parents=True, exist_ok=True)
                
                # Execute actual moves
                for operation in operations:
                    move_result = self._execute_single_move(operation)
                    
                    if move_result.success:
                        self.move_stats["directories_moved"] += 1
                        self.move_stats["total_size_moved_mb"] += operation.size_bytes / (1024 * 1024)
                        
                        if operation.backup_created:
                            self.move_stats["backups_created"] += 1
                    else:
                        self.move_stats["directories_failed"] += 1
                        result.errors.extend(move_result.errors)
            
            # Add operation details to result
            result.metadata.update({
                "operations": [
                    {
                        "source": op.source_name,
                        "destination": str(op.destination_path),
                        "status": op.status,
                        "size_mb": op.size_bytes / (1024 * 1024),
                        "error": op.error_message
                    }
                    for op in operations
                ]
            })
            
        except Exception as e:
            result.add_error(f"Failed to execute move operations: {str(e)}")
        
        return result
    
    def _execute_single_move(self, operation: MoveOperation) -> ProcessingResult:
        """Execute single move operation with retries."""
        result = ProcessingResult(success=True)
        
        for attempt in range(self.move_config.max_retries):
            try:
                self.logger.debug(f"Attempting move", 
                                source=str(operation.source_path),
                                destination=str(operation.destination_path),
                                attempt=attempt + 1)
                
                # Create backup if configured
                if self.move_config.create_backup and operation.destination_path.exists():
                    backup_path = operation.destination_path.with_suffix(
                        operation.destination_path.suffix + self.move_config.backup_suffix
                    )
                    shutil.move(str(operation.destination_path), str(backup_path))
                    operation.backup_created = True
                
                # Perform move
                shutil.move(str(operation.source_path), str(operation.destination_path))
                
                # Verify move if configured
                if self.move_config.verify_moves:
                    if not operation.destination_path.exists():
                        raise FileNotFoundError(f"Move verification failed: {operation.destination_path}")
                    if operation.source_path.exists():
                        raise FileExistsError(f"Source still exists after move: {operation.source_path}")
                
                # Preserve permissions if configured
                if self.move_config.preserve_permissions:
                    try:
                        # Copy permissions from source parent (best effort)
                        source_stat = operation.source_path.parent.stat()
                        os.chmod(operation.destination_path, source_stat.st_mode)
                    except Exception:
                        # Ignore permission errors
                        pass
                
                operation.status = "success"
                
                self.logger.info(f"Directory moved successfully", 
                               source=operation.source_name,
                               destination=str(operation.destination_path))
                
                break  # Success, exit retry loop
                
            except Exception as e:
                operation.error_message = str(e)
                
                if attempt == self.move_config.max_retries - 1:
                    # Final attempt failed
                    operation.status = "failed"
                    result.add_error(f"Failed to move {operation.source_name}: {str(e)}")
                    
                    self.logger.error(f"Failed to move directory", 
                                    source=operation.source_name,
                                    error=str(e),
                                    attempts=attempt + 1)
                else:
                    # Retry
                    self.logger.warning(f"Move attempt failed, retrying", 
                                      source=operation.source_name,
                                      error=str(e),
                                      attempt=attempt + 1)
        
        return result
    
    def get_move_statistics(self) -> Dict[str, Any]:
        """Get comprehensive move statistics."""
        base_stats = self.get_statistics()
        
        success_rate = 0.0
        if self.move_stats["directories_found"] > 0:
            success_rate = self.move_stats["directories_moved"] / self.move_stats["directories_found"]
        
        return {
            **base_stats,
            **self.move_stats,
            "success_rate": success_rate
        }
    
    def validate_input(self, input_data: Path) -> bool:
        """Validate input directory."""
        try:
            input_path = Path(input_data) if isinstance(input_data, str) else input_data
            return input_path.exists() and input_path.is_dir()
        except Exception:
            return False


def main():
    """Main entry point for directory moving."""
    # Configuration
    source_dir = Path("/groups/itay_mayrose/alongonda/datasets/full_genomes/plaza")
    destination_dir = Path("/groups/itay_mayrose/alongonda/datasets/full_genomes/plaza/organisms")
    exclude_dirs = ["processed_annotations"]
    
    # Setup configuration
    config = BioinformaticsConfig(
        input_dir=source_dir,
        output_dir=destination_dir,
        max_workers=1,
        log_level="INFO"
    )
    
    move_config = DirectoryMoveConfig(
        exclude_patterns=[],
        create_backup=False,
        dry_run=False,
        overwrite_existing=False,
        verify_moves=True
    )
    
    # Initialize mover
    try:
        mover = DirectoryMover(
            config=config,
            move_config=move_config
        )
        
        # Run moving
        result = mover.run(source_dir, 
                          destination_dir=destination_dir,
                          exclude_directories=exclude_dirs)
        
        # Display results
        if result.success:
            print("✅ Directory moving completed successfully!")
            stats = mover.get_move_statistics()
            print(f"   Directories found: {stats['directories_found']}")
            print(f"   Directories moved: {stats['directories_moved']}")
            print(f"   Directories failed: {stats['directories_failed']}")
            print(f"   Total size moved: {stats['total_size_moved_mb']:.2f} MB")
            print(f"   Success rate: {stats['success_rate']:.2%}")
            print(f"   Processing time: {result.processing_time:.2f} seconds")
            
            if move_config.dry_run:
                print("\n🔍 DRY RUN MODE - No actual changes were made")
        else:
            print("❌ Directory moving failed!")
            for error in result.errors:
                print(f"   Error: {error}")
        
        # Display warnings
        for warning in result.warnings:
            print(f"   Warning: {warning}")
        
    except Exception as e:
        print(f"❌ Critical error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())