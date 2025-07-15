"""
Industry-Level File Renamer with Species Mapping.

This module provides object-oriented file renaming functionality with comprehensive
error handling, validation, backup capabilities, and batch processing.
"""

from ..config.settings import get_settings
from ..core.base import BioinformaticsProcessor
from ..core.exceptions import ProcessingError, ValidationError
from ..utils.file_operations import file_manager
from loguru import logger
import os
import shutil
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from ..base_classes import (
    BaseProcessor,
    BioinformaticsConfig,
    ProcessingResult,
    FileSystemProcessor
)


@dataclass
class RenamingRule:
    """Configuration for a single renaming rule."""
    
    old_pattern: str
    new_pattern: str
    case_sensitive: bool = True
    exact_match: bool = True
    
    def apply_transformation(self, text: str) -> str:
        """Apply transformation to text."""
        if not self.case_sensitive:
            text = text.lower()
        
        if self.exact_match:
            return self.new_pattern if text == self.old_pattern else text
        else:
            return text.replace(self.old_pattern, self.new_pattern)


@dataclass
class FileRenamingConfig:
    """Configuration for file renaming operations."""
    
    file_extensions: List[str] = field(default_factory=lambda: [".csv"])
    preserve_extensions: bool = True
    create_backup: bool = False
    backup_suffix: str = ".backup"
    dry_run: bool = False
    case_sensitive: bool = False
    normalize_names: bool = True
    replace_spaces: bool = True
    space_replacement: str = "_"
    remove_special_chars: bool = False
    max_filename_length: int = 255
    
    def normalize_filename(self, filename: str) -> str:
        """Normalize filename according to configuration."""
        if self.normalize_names:
            filename = filename.lower()
        
        if self.replace_spaces:
            filename = filename.replace(" ", self.space_replacement)
        
        if self.remove_special_chars:
            # Keep only alphanumeric, underscore, hyphen, and dot
            filename = "".join(c for c in filename if c.isalnum() or c in "._-")
        
        # Truncate if too long
        if len(filename) > self.max_filename_length:
            name, ext = os.path.splitext(filename)
            max_name_length = self.max_filename_length - len(ext)
            filename = name[:max_name_length] + ext
        
        return filename


@dataclass
class RenameOperation:
    """Information about a file rename operation."""
    
    original_path: Path
    new_path: Path
    rule_applied: Optional[str] = None
    backup_created: bool = False
    status: str = "pending"  # pending, success, failed, skipped
    error_message: Optional[str] = None
    
    @property
    def original_name(self) -> str:
        """Get original filename."""
        return self.original_path.name
    
    @property
    def new_name(self) -> str:
        """Get new filename."""
        return self.new_path.name


class SpeciesMapLoader:
    """Loader for species mapping data."""
    
    def __init__(self, logger):
        """Initialize species map loader."""
        self.logger = logger
    
    def load_species_mapping(self, mapping_file: Path, delimiter: str = '\t') -> Dict[str, str]:
        """
        Load species mapping from file.
        
        Args:
            mapping_file: Path to mapping file
            delimiter: File delimiter
            
        Returns:
            Dictionary mapping old names to new names
        """
        try:
            self.logger.info(f"Loading species mapping", file=str(mapping_file))
            
            df = pd.read_csv(mapping_file, delimiter=delimiter)
            
            if len(df.columns) < 2:
                raise ValueError(f"Mapping file must have at least 2 columns, got {len(df.columns)}")
            
            # Use first two columns as key-value pairs
            keys = df.iloc[:, 0].astype(str).values
            values = df.iloc[:, 1].astype(str).values
            
            species_map = dict(zip(keys, values))
            
            self.logger.info(f"Loaded species mapping", 
                           entries=len(species_map),
                           sample_entries=dict(list(species_map.items())[:3]))
            
            return species_map
            
        except Exception as e:
            self.logger.error(f"Failed to load species mapping", 
                            file=str(mapping_file),
                            error=str(e))
            raise


class FileRenamer(FileSystemProcessor):
    """
    Professional file renamer with comprehensive functionality.
    """
    
    def __init__(self, 
                 config: BioinformaticsConfig,
                 renaming_config: FileRenamingConfig):
        """
        Initialize file renamer.
        
        Args:
            config: Base processing configuration
            renaming_config: Renaming-specific configuration
        """
        super().__init__(config, name="FileRenamer")
        
        self.renaming_config = renaming_config
        self.species_loader = SpeciesMapLoader(self.logger)
        
        # Statistics
        self.rename_stats = {
            "files_found": 0,
            "files_renamed": 0,
            "files_skipped": 0,
            "files_failed": 0,
            "backups_created": 0,
            "rules_applied": 0,
        }
    
    def process(self, input_data: Path, **kwargs) -> ProcessingResult:
        """
        Process file renaming operations.
        
        Args:
            input_data: Path to directory containing files to rename
            **kwargs: Additional processing parameters including:
                - species_mapping_file: Path to species mapping file
                - custom_mapping: Dictionary of custom name mappings
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            input_dir = Path(input_data) if isinstance(input_data, str) else input_data
            
            self.logger.info(f"Starting file renaming", 
                           input_dir=str(input_dir),
                           dry_run=self.renaming_config.dry_run)
            
            # Load mappings
            mappings = {}
            
            if 'species_mapping_file' in kwargs:
                species_mapping = self.species_loader.load_species_mapping(
                    Path(kwargs['species_mapping_file'])
                )
                mappings.update(species_mapping)
            
            if 'custom_mapping' in kwargs:
                mappings.update(kwargs['custom_mapping'])
            
            if not mappings:
                result.add_error("No renaming mappings provided")
                return result
            
            # Step 1: Find files to rename
            rename_operations = self._plan_rename_operations(input_dir, mappings)
            self.rename_stats["files_found"] = len(rename_operations)
            
            if not rename_operations:
                result.add_warning("No files found matching renaming criteria")
                return result
            
            # Step 2: Execute renaming operations
            execution_result = self._execute_rename_operations(rename_operations)
            result = result.merge(execution_result)
            
            # Update metadata
            result.metadata.update({
                "rename_stats": self.rename_stats,
                "input_directory": str(input_dir),
                "dry_run": self.renaming_config.dry_run,
                "mappings_count": len(mappings),
                "operations_planned": len(rename_operations)
            })
            
            self.logger.info(f"File renaming completed",
                           files_renamed=self.rename_stats["files_renamed"],
                           files_failed=self.rename_stats["files_failed"])
            
        except Exception as e:
            result.add_error(f"File renaming failed: {str(e)}")
        
        return result
    
    def _plan_rename_operations(self, input_dir: Path, mappings: Dict[str, str]) -> List[RenameOperation]:
        """Plan rename operations for files in directory."""
        operations = []
        
        try:
            for file_path in input_dir.iterdir():
                if not file_path.is_file():
                    continue
                
                # Check if file extension matches
                if not any(file_path.name.lower().endswith(ext.lower()) 
                          for ext in self.renaming_config.file_extensions):
                    continue
                
                # Get base name without extension
                base_name = file_path.stem
                extension = file_path.suffix if self.renaming_config.preserve_extensions else ""
                
                # Check for mapping
                new_base_name = None
                rule_applied = None
                
                # Try exact match first
                if base_name in mappings:
                    new_base_name = mappings[base_name]
                    rule_applied = f"exact_match: {base_name} -> {new_base_name}"
                elif not self.renaming_config.case_sensitive:
                    # Try case-insensitive match
                    for old_name, new_name in mappings.items():
                        if base_name.lower() == old_name.lower():
                            new_base_name = new_name
                            rule_applied = f"case_insensitive_match: {old_name} -> {new_name}"
                            break
                
                if new_base_name:
                    # Normalize new filename
                    normalized_name = self.renaming_config.normalize_filename(new_base_name)
                    new_filename = f"{normalized_name}{extension}"
                    new_path = input_dir / new_filename
                    
                    # Skip if new path already exists and is different from current
                    if new_path.exists() and new_path != file_path:
                        self.logger.warning(f"Target file already exists", 
                                          original=str(file_path),
                                          target=str(new_path))
                        continue
                    
                    # Skip if no actual change needed
                    if new_path == file_path:
                        self.logger.debug(f"No change needed", file=str(file_path))
                        continue
                    
                    operation = RenameOperation(
                        original_path=file_path,
                        new_path=new_path,
                        rule_applied=rule_applied
                    )
                    operations.append(operation)
                    
                    self.logger.debug(f"Planned rename operation", 
                                    original=file_path.name,
                                    new=new_filename,
                                    rule=rule_applied)
            
            self.logger.info(f"Planned rename operations", count=len(operations))
            return operations
            
        except Exception as e:
            self.logger.error(f"Failed to plan rename operations", error=str(e))
            raise
    
    def _execute_rename_operations(self, operations: List[RenameOperation]) -> ProcessingResult:
        """Execute planned rename operations."""
        result = ProcessingResult(success=True)
        
        try:
            if self.renaming_config.dry_run:
                self.logger.info(f"DRY RUN: Would execute rename operations", count=len(operations))
                
                for operation in operations:
                    operation.status = "dry_run"
                    self.rename_stats["files_renamed"] += 1
                    
                    self.logger.info(f"DRY RUN: Would rename", 
                                   original=operation.original_name,
                                   new=operation.new_name)
            else:
                # Execute actual renaming
                for operation in operations:
                    try:
                        # Create backup if configured
                        if self.renaming_config.create_backup:
                            backup_path = operation.original_path.with_suffix(
                                operation.original_path.suffix + self.renaming_config.backup_suffix
                            )
                            shutil.copy2(operation.original_path, backup_path)
                            operation.backup_created = True
                            self.rename_stats["backups_created"] += 1
                        
                        # Perform rename
                        operation.original_path.rename(operation.new_path)
                        operation.status = "success"
                        self.rename_stats["files_renamed"] += 1
                        self.rename_stats["rules_applied"] += 1
                        
                        self.logger.info(f"File renamed successfully", 
                                       original=operation.original_name,
                                       new=operation.new_name)
                        
                    except Exception as e:
                        operation.status = "failed"
                        operation.error_message = str(e)
                        self.rename_stats["files_failed"] += 1
                        
                        result.add_error(f"Failed to rename {operation.original_name}: {str(e)}")
                        self.logger.error(f"Failed to rename file", 
                                        file=operation.original_name,
                                        error=str(e))
            
            # Add operation details to result
            result.metadata.update({
                "operations": [
                    {
                        "original": op.original_name,
                        "new": op.new_name,
                        "status": op.status,
                        "rule": op.rule_applied,
                        "error": op.error_message
                    }
                    for op in operations
                ]
            })
            
        except Exception as e:
            result.add_error(f"Failed to execute rename operations: {str(e)}")
        
        return result
    
    def get_rename_statistics(self) -> Dict[str, Any]:
        """Get comprehensive rename statistics."""
        base_stats = self.get_statistics()
        
        success_rate = 0.0
        if self.rename_stats["files_found"] > 0:
            success_rate = self.rename_stats["files_renamed"] / self.rename_stats["files_found"]
        
        return {
            **base_stats,
            **self.rename_stats,
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
    """Main entry point for file renaming."""
    # Configuration
    input_dir = Path("/groups/itay_mayrose/alongonda/datasets/full_genomes/plaza/processed_annotations")
    species_mapping_file = Path("/groups/itay_mayrose/alongonda/desktop/species_information.csv")
    
    # Setup configuration
    config = BioinformaticsConfig(
        input_dir=input_dir,
        output_dir=input_dir,
        max_workers=1,
        log_level="INFO"
    )
    
    renaming_config = FileRenamingConfig(
        file_extensions=[".csv"],
        preserve_extensions=True,
        create_backup=False,
        dry_run=False,
        normalize_names=True,
        replace_spaces=True,
        space_replacement="_"
    )
    
    # Initialize renamer
    try:
        renamer = FileRenamer(
            config=config,
            renaming_config=renaming_config
        )
        
        # Run renaming
        result = renamer.run(input_dir, species_mapping_file=species_mapping_file)
        
        # Display results
        if result.success:
            logger.info("✅ File renaming completed successfully!")
            stats = renamer.get_rename_statistics()
            logger.info(f"   Files found: {stats['files_found']}")
            logger.info(f"   Files renamed: {stats['files_renamed']}")
            logger.info(f"   Files failed: {stats['files_failed']}")
            logger.info(f"   Success rate: {stats['success_rate']:.2%}")
            logger.info(f"   Processing time: {result.processing_time:.2f} seconds")
            
            if renaming_config.dry_run:
                logger.info("\n🔍 DRY RUN MODE - No actual changes were made")
        else:
            logger.info("❌ File renaming failed!")
            for error in result.errors:
                logger.info(f"   Error: {error}")
        
        # Display warnings
        for warning in result.warnings:
            logger.info(f"   Warning: {warning}")
        
    except Exception as e:
        logger.info(f"❌ Critical error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())