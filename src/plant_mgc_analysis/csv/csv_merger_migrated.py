"""
Industry-Level CSV File Merger.

This module provides object-oriented CSV file merging with comprehensive
error handling, validation, and parallel processing capabilities.
"""

from ..config.settings import get_settings
from ..core.base import BioinformaticsProcessor
from ..core.exceptions import ProcessingError, ValidationError
from ..utils.file_operations import file_manager
from loguru import logger
import os
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from ..base_classes import (
    BaseProcessor,
    BioinformaticsConfig,
    ProcessingResult,
    DatabaseProcessor,
    FileSystemProcessor
)


@dataclass
class CSVMergeConfig:
    """Configuration for CSV merging operations."""
    
    pathway_name_extraction: str = "filename_prefix"  # "filename_prefix", "directory_name", "custom"
    separator_char: str = "_"
    pathway_column_name: str = "pathway"
    validate_schemas: bool = True
    skip_empty_files: bool = True
    max_memory_usage_mb: int = 1000
    chunk_size: int = 10000
    preserve_dtypes: bool = True
    
    def extract_pathway_name(self, file_path: Path) -> str:
        """Extract pathway name based on configuration."""
        if self.pathway_name_extraction == "filename_prefix":
            return file_path.stem.split(self.separator_char)[0]
        elif self.pathway_name_extraction == "directory_name":
            return file_path.parent.name
        else:
            return file_path.stem


@dataclass
class CSVFileInfo:
    """Information about a CSV file to be merged."""
    
    file_path: Path
    pathway_name: str
    row_count: int = 0
    column_count: int = 0
    file_size_mb: float = 0.0
    columns: List[str] = None
    dtypes: Dict[str, Any] = None
    
    def __post_init__(self):
        """Calculate file statistics."""
        if self.file_path.exists():
            self.file_size_mb = self.file_path.stat().st_size / (1024 * 1024)


class CSVSchemaValidator:
    """Validator for CSV schema consistency."""
    
    def __init__(self, config: CSVMergeConfig, logger):
        """Initialize schema validator."""
        self.config = config
        self.logger = logger
        self.reference_schema = None
    
    def validate_file_schema(self, file_info: CSVFileInfo) -> Tuple[bool, List[str]]:
        """
        Validate CSV file schema consistency.
        
        Args:
            file_info: CSV file information
            
        Returns:
            Tuple of (is_valid, error_messages)
        """
        errors = []
        
        try:
            # Read first few rows to check schema
            df_sample = pd.read_csv(file_info.file_path, nrows=5)
            file_info.columns = list(df_sample.columns)
            file_info.column_count = len(df_sample.columns)
            file_info.dtypes = df_sample.dtypes.to_dict()
            
            if self.reference_schema is None:
                # First file - set as reference
                self.reference_schema = {
                    'columns': set(file_info.columns),
                    'column_count': file_info.column_count
                }
                self.logger.info(f"Set reference schema", 
                               file=str(file_info.file_path),
                               columns=file_info.column_count)
            else:
                # Validate against reference
                if set(file_info.columns) != self.reference_schema['columns']:
                    errors.append(f"Column mismatch in {file_info.file_path}: "
                                f"expected {self.reference_schema['columns']}, "
                                f"got {set(file_info.columns)}")
                
                if file_info.column_count != self.reference_schema['column_count']:
                    errors.append(f"Column count mismatch in {file_info.file_path}: "
                                f"expected {self.reference_schema['column_count']}, "
                                f"got {file_info.column_count}")
            
            return len(errors) == 0, errors
            
        except Exception as e:
            errors.append(f"Failed to validate schema for {file_info.file_path}: {str(e)}")
            return False, errors


class CSVFileMerger(FileSystemProcessor):
    """
    Professional CSV file merger with comprehensive functionality.
    """
    
    def __init__(self, 
                 config: BioinformaticsConfig,
                 merge_config: CSVMergeConfig):
        """
        Initialize CSV file merger.
        
        Args:
            config: Base processing configuration
            merge_config: Merge-specific configuration
        """
        super().__init__(config, name="CSVFileMerger")
        
        self.merge_config = merge_config
        self.schema_validator = CSVSchemaValidator(merge_config, self.logger)
        
        # Statistics
        self.merge_stats = {
            "files_found": 0,
            "files_processed": 0,
            "files_failed": 0,
            "total_rows_merged": 0,
            "total_size_mb": 0.0,
        }
    
    def process(self, input_data: Path, **kwargs) -> ProcessingResult:
        """
        Process CSV file merging.
        
        Args:
            input_data: Path to directory containing CSV files
            **kwargs: Additional processing parameters including:
                - output_file: Output CSV file path
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            input_dir = Path(input_data) if isinstance(input_data, str) else input_data
            output_file = kwargs.get('output_file', input_dir / "merged_pathways.csv")
            
            self.logger.info(f"Starting CSV file merging", 
                           input_dir=str(input_dir),
                           output_file=str(output_file))
            
            # Step 1: Find and analyze CSV files
            csv_files = self._find_csv_files(input_dir)
            self.merge_stats["files_found"] = len(csv_files)
            
            if not csv_files:
                result.add_error("No CSV files found in directory")
                return result
            
            # Step 2: Validate schemas if configured
            if self.merge_config.validate_schemas:
                validation_result = self._validate_schemas(csv_files)
                if not validation_result.success:
                    result.errors.extend(validation_result.errors)
                    return result
            
            # Step 3: Merge files
            merge_result = self._merge_csv_files(csv_files, Path(output_file))
            result = result.merge(merge_result)
            
            # Update metadata
            result.metadata.update({
                "merge_stats": self.merge_stats,
                "input_directory": str(input_dir),
                "output_file": str(output_file),
                "files_merged": len([f for f in csv_files if f.row_count > 0])
            })
            
            self.logger.info(f"CSV merging completed",
                           files_processed=self.merge_stats["files_processed"],
                           total_rows=self.merge_stats["total_rows_merged"])
            
        except Exception as e:
            result.add_error(f"CSV merging failed: {str(e)}")
        
        return result
    
    def _find_csv_files(self, input_dir: Path) -> List[CSVFileInfo]:
        """Find and analyze CSV files in directory."""
        csv_files = []
        
        try:
            for root, dirs, files in os.walk(input_dir):
                for file_name in files:
                    if file_name.endswith('.csv'):
                        file_path = Path(root) / file_name
                        
                        if not file_path.is_file():
                            continue
                        
                        # Create file info
                        pathway_name = self.merge_config.extract_pathway_name(file_path)
                        file_info = CSVFileInfo(
                            file_path=file_path,
                            pathway_name=pathway_name
                        )
                        
                        # Skip empty files if configured
                        if self.merge_config.skip_empty_files and file_info.file_size_mb < 0.001:
                            self.logger.warning(f"Skipping empty file", file=str(file_path))
                            continue
                        
                        csv_files.append(file_info)
                        self.merge_stats["total_size_mb"] += file_info.file_size_mb
                        
                        self.logger.debug(f"Found CSV file", 
                                        file=str(file_path),
                                        pathway=pathway_name,
                                        size_mb=file_info.file_size_mb)
            
            self.logger.info(f"Found CSV files", 
                           count=len(csv_files),
                           total_size_mb=self.merge_stats["total_size_mb"])
            
            return csv_files
            
        except Exception as e:
            self.logger.error(f"Failed to find CSV files", error=str(e))
            raise
    
    def _validate_schemas(self, csv_files: List[CSVFileInfo]) -> ProcessingResult:
        """Validate CSV file schemas for consistency."""
        result = ProcessingResult(success=True)
        
        try:
            self.logger.info(f"Validating CSV schemas", file_count=len(csv_files))
            
            for file_info in csv_files:
                is_valid, errors = self.schema_validator.validate_file_schema(file_info)
                
                if not is_valid:
                    for error in errors:
                        result.add_error(error)
                        self.logger.error(error)
                else:
                    self.logger.debug(f"Schema validation passed", 
                                    file=str(file_info.file_path))
            
            return result
            
        except Exception as e:
            result.add_error(f"Schema validation failed: {str(e)}")
            return result
    
    def _merge_csv_files(self, csv_files: List[CSVFileInfo], output_file: Path) -> ProcessingResult:
        """Merge CSV files into single output file."""
        result = ProcessingResult(success=True)
        
        try:
            dataframes = []
            
            # Process files in parallel or sequential based on configuration
            if self.config.max_workers == 1:
                # Sequential processing
                for file_info in csv_files:
                    df_result = self._process_single_csv(file_info)
                    if df_result:
                        dataframes.append(df_result)
            else:
                # Parallel processing
                with ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
                    future_to_file = {
                        executor.submit(self._process_single_csv, file_info): file_info
                        for file_info in csv_files
                    }
                    
                    for future in as_completed(future_to_file):
                        file_info = future_to_file[future]
                        try:
                            df_result = future.result(timeout=self.config.timeout)
                            if df_result is not None:
                                dataframes.append(df_result)
                        except Exception as e:
                            result.add_error(f"Failed to process {file_info.file_path}: {str(e)}")
                            self.merge_stats["files_failed"] += 1
            
            if not dataframes:
                result.add_error("No valid dataframes to merge")
                return result
            
            # Merge all dataframes
            self.logger.info(f"Merging dataframes", count=len(dataframes))
            merged_df = pd.concat(dataframes, ignore_index=True)
            
            # Save merged result
            output_file.parent.mkdir(parents=True, exist_ok=True)
            merged_df.to_csv(output_file, index=False)
            
            # Update statistics
            self.merge_stats["total_rows_merged"] = len(merged_df)
            
            result.metadata.update({
                "output_file": str(output_file),
                "total_rows": len(merged_df),
                "total_columns": len(merged_df.columns),
                "dataframes_merged": len(dataframes)
            })
            
            self.logger.info(f"CSV files merged successfully", 
                           output=str(output_file),
                           total_rows=len(merged_df),
                           dataframes=len(dataframes))
            
        except Exception as e:
            result.add_error(f"Failed to merge CSV files: {str(e)}")
        
        return result
    
    def _process_single_csv(self, file_info: CSVFileInfo) -> Optional[pd.DataFrame]:
        """Process single CSV file and add pathway information."""
        try:
            self.logger.debug(f"Processing CSV file", file=str(file_info.file_path))
            
            # Read CSV file
            df = pd.read_csv(file_info.file_path)
            
            # Add pathway column
            df[self.merge_config.pathway_column_name] = file_info.pathway_name
            
            # Update file statistics
            file_info.row_count = len(df)
            
            # Update global statistics
            self.merge_stats["files_processed"] += 1
            
            self.logger.debug(f"Processed CSV file successfully", 
                            file=str(file_info.file_path),
                            rows=len(df),
                            pathway=file_info.pathway_name)
            
            return df
            
        except Exception as e:
            self.logger.error(f"Failed to process CSV file", 
                            file=str(file_info.file_path),
                            error=str(e))
            self.merge_stats["files_failed"] += 1
            return None
    
    def get_merge_statistics(self) -> Dict[str, Any]:
        """Get comprehensive merge statistics."""
        base_stats = self.get_statistics()
        
        success_rate = 0.0
        if self.merge_stats["files_found"] > 0:
            success_rate = self.merge_stats["files_processed"] / self.merge_stats["files_found"]
        
        return {
            **base_stats,
            **self.merge_stats,
            "success_rate": success_rate,
            "average_rows_per_file": (
                self.merge_stats["total_rows_merged"] / 
                max(1, self.merge_stats["files_processed"])
            )
        }
    
    def validate_input(self, input_data: Path) -> bool:
        """Validate input directory."""
        try:
            input_path = Path(input_data) if isinstance(input_data, str) else input_data
            return input_path.exists() and input_path.is_dir()
        except Exception:
            return False


def main():
    """Main entry point for CSV file merging."""
    # Configuration
    input_dir = Path("/groups/itay_mayrose/alongonda/datasets/KEGG_fasta_updated_fixed")
    output_file = input_dir / "merged_pathways.csv"
    
    # Setup configuration
    config = BioinformaticsConfig(
        input_dir=input_dir,
        output_dir=input_dir,
        max_workers=4,
        log_level="INFO"
    )
    
    merge_config = CSVMergeConfig(
        pathway_name_extraction="filename_prefix",
        separator_char="_",
        validate_schemas=True,
        skip_empty_files=True
    )
    
    # Initialize merger
    try:
        merger = CSVFileMerger(
            config=config,
            merge_config=merge_config
        )
        
        # Run merging
        result = merger.run(input_dir, output_file=output_file)
        
        # Display results
        if result.success:
            logger.info("✅ CSV file merging completed successfully!")
            stats = merger.get_merge_statistics()
            logger.info(f"   Files processed: {stats['files_processed']}")
            logger.info(f"   Total rows merged: {stats['total_rows_merged']}")
            logger.info(f"   Success rate: {stats['success_rate']:.2%}")
            logger.info(f"   Processing time: {result.processing_time:.2f} seconds")
            logger.info(f"   Output file: {result.metadata.get('output_file')}")
        else:
            logger.info("❌ CSV file merging failed!")
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
