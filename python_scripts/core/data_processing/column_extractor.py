"""
Industry-Level Column Extractor and Data Processor.

This module provides object-oriented CSV column extraction and processing with comprehensive
error handling, validation, data cleaning, and parallel processing capabilities.
"""

import csv
import re
import sys
import os
from typing import Dict, List, Optional, Any, Tuple, Callable
from dataclasses import dataclass, field
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from ..base_classes import (
    BaseProcessor,
    BioinformaticsConfig,
    ProcessingResult,
    DatabaseProcessor,
    ParallelProcessor
)


@dataclass
class ColumnMapping:
    """Configuration for column mapping and extraction."""
    
    source_columns: List[str] = field(default_factory=lambda: ["id", "", "", "start", "end", "sequence"])
    target_columns: List[str] = field(default_factory=lambda: ["id", "start", "end", "sequence"])
    column_indices: List[int] = field(default_factory=lambda: [0, 3, 4, 5])
    required_columns: int = 6
    header_row: bool = True


@dataclass
class DataCleaningConfig:
    """Configuration for data cleaning operations."""
    
    remove_parentheses: bool = True
    parentheses_pattern: str = r'\(.*?\)'
    normalize_whitespace: bool = True
    remove_empty_rows: bool = True
    validate_data_types: bool = True
    max_field_size: Optional[int] = None
    
    def get_cleaning_functions(self) -> List[Callable[[str], str]]:
        """Get list of cleaning functions to apply."""
        functions = []
        
        if self.remove_parentheses:
            functions.append(lambda text: re.sub(self.parentheses_pattern, '', text))
        
        if self.normalize_whitespace:
            functions.append(lambda text: ' '.join(text.split()))
        
        return functions


@dataclass
class ProcessingStats:
    """Statistics for processing operations."""
    
    files_processed: int = 0
    files_failed: int = 0
    rows_processed: int = 0
    rows_extracted: int = 0
    rows_failed: int = 0
    data_cleaning_applied: int = 0


class CSVFieldSizeManager:
    """Manager for CSV field size limits."""
    
    def __init__(self, logger):
        """Initialize CSV field size manager."""
        self.logger = logger
        self.original_limit = csv.field_size_limit()
    
    def increase_field_size_limit(self, max_size: Optional[int] = None) -> int:
        """Increase CSV field size limit safely."""
        try:
            target_size = max_size or sys.maxsize
            
            while True:
                try:
                    csv.field_size_limit(target_size)
                    self.logger.debug(f"Set CSV field size limit", size=target_size)
                    return target_size
                except OverflowError:
                    target_size = int(target_size / 10)
                    if target_size < 1000:
                        # Fallback to a reasonable minimum
                        csv.field_size_limit(131072)  # 128KB
                        return 131072
        except Exception as e:
            self.logger.warning(f"Failed to set field size limit", error=str(e))
            return self.original_limit
    
    def restore_original_limit(self):
        """Restore original CSV field size limit."""
        try:
            csv.field_size_limit(self.original_limit)
        except Exception:
            pass


class DataCleaner:
    """Data cleaning and validation utilities."""
    
    def __init__(self, config: DataCleaningConfig, logger):
        """Initialize data cleaner."""
        self.config = config
        self.logger = logger
        self.cleaning_functions = config.get_cleaning_functions()
    
    def clean_text(self, text: str) -> str:
        """Apply all configured cleaning functions to text."""
        if not text:
            return text
        
        cleaned_text = text
        for func in self.cleaning_functions:
            try:
                cleaned_text = func(cleaned_text)
            except Exception as e:
                self.logger.warning(f"Text cleaning function failed", 
                                  text=text[:50],
                                  error=str(e))
        
        return cleaned_text
    
    def validate_row_data(self, row_data: List[str], expected_columns: int) -> Tuple[bool, List[str]]:
        """Validate extracted row data."""
        errors = []
        
        if len(row_data) < expected_columns:
            errors.append(f"Insufficient columns: expected {expected_columns}, got {len(row_data)}")
        
        # Check for empty required fields
        for i, value in enumerate(row_data):
            if not value or value.strip() == "":
                errors.append(f"Empty value in column {i}")
        
        # Validate start/end positions if they exist
        if len(row_data) >= 3:
            try:
                start_pos = int(row_data[1]) if row_data[1].isdigit() else None
                end_pos = int(row_data[2]) if row_data[2].isdigit() else None
                
                if start_pos is not None and end_pos is not None and start_pos >= end_pos:
                    errors.append(f"Invalid coordinates: start ({start_pos}) >= end ({end_pos})")
            except (ValueError, IndexError):
                pass  # Non-numeric coordinates are acceptable in some cases
        
        return len(errors) == 0, errors


class ColumnExtractor(DatabaseProcessor):
    """
    Professional column extractor with comprehensive functionality.
    """
    
    def __init__(self, 
                 config: BioinformaticsConfig,
                 column_mapping: ColumnMapping,
                 cleaning_config: DataCleaningConfig):
        """
        Initialize column extractor.
        
        Args:
            config: Base processing configuration
            column_mapping: Column mapping configuration
            cleaning_config: Data cleaning configuration
        """
        super().__init__(config, name="ColumnExtractor")
        
        self.column_mapping = column_mapping
        self.cleaning_config = cleaning_config
        self.field_manager = CSVFieldSizeManager(self.logger)
        self.data_cleaner = DataCleaner(cleaning_config, self.logger)
        
        # Statistics
        self.stats = ProcessingStats()
    
    def process(self, input_data: Path, **kwargs) -> ProcessingResult:
        """
        Process column extraction from CSV files.
        
        Args:
            input_data: Path to input directory or single CSV file
            **kwargs: Additional processing parameters including:
                - output_dir: Output directory path
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            input_path = Path(input_data) if isinstance(input_data, str) else input_data
            output_dir = Path(kwargs.get('output_dir', self.config.output_dir))
            
            self.logger.info(f"Starting column extraction", 
                           input=str(input_path),
                           output=str(output_dir))
            
            # Setup CSV field size limit
            self.field_manager.increase_field_size_limit(self.cleaning_config.max_field_size)
            
            try:
                # Determine if input is file or directory
                if input_path.is_file():
                    if input_path.suffix.lower() == '.csv':
                        extraction_result = self._process_single_file(input_path, output_dir)
                        result = result.merge(extraction_result)
                    else:
                        result.add_error(f"Input file is not a CSV: {input_path}")
                
                elif input_path.is_dir():
                    directory_result = self._process_directory(input_path, output_dir)
                    result = result.merge(directory_result)
                
                else:
                    result.add_error(f"Input path does not exist: {input_path}")
                
            finally:
                # Restore original CSV field size limit
                self.field_manager.restore_original_limit()
            
            # Update metadata
            result.metadata.update({
                "processing_stats": {
                    "files_processed": self.stats.files_processed,
                    "files_failed": self.stats.files_failed,
                    "rows_processed": self.stats.rows_processed,
                    "rows_extracted": self.stats.rows_extracted,
                    "rows_failed": self.stats.rows_failed,
                    "data_cleaning_applied": self.stats.data_cleaning_applied
                },
                "input_path": str(input_path),
                "output_directory": str(output_dir)
            })
            
            self.logger.info(f"Column extraction completed",
                           files_processed=self.stats.files_processed,
                           rows_extracted=self.stats.rows_extracted)
            
        except Exception as e:
            result.add_error(f"Column extraction failed: {str(e)}")
        
        return result
    
    def _process_directory(self, input_dir: Path, output_dir: Path) -> ProcessingResult:
        """Process all CSV files in directory."""
        result = ProcessingResult(success=True)
        
        try:
            # Find CSV files
            csv_files = []
            for root, dirs, files in os.walk(input_dir):
                for filename in files:
                    if filename.endswith(".csv"):
                        csv_files.append(Path(root) / filename)
            
            if not csv_files:
                result.add_warning("No CSV files found in directory")
                return result
            
            self.logger.info(f"Found CSV files", count=len(csv_files))
            
            # Process files in parallel or sequential
            if self.config.max_workers == 1:
                # Sequential processing
                for csv_file in csv_files:
                    file_result = self._process_single_file(csv_file, output_dir, input_dir)
                    if not file_result.success:
                        result.errors.extend(file_result.errors)
                        result.warnings.extend(file_result.warnings)
            else:
                # Parallel processing
                with ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
                    future_to_file = {
                        executor.submit(self._process_single_file, csv_file, output_dir, input_dir): csv_file
                        for csv_file in csv_files
                    }
                    
                    for future in as_completed(future_to_file):
                        csv_file = future_to_file[future]
                        try:
                            file_result = future.result(timeout=self.config.timeout)
                            if not file_result.success:
                                result.errors.extend(file_result.errors)
                                result.warnings.extend(file_result.warnings)
                        except Exception as e:
                            result.add_error(f"Failed to process {csv_file}: {str(e)}")
                            self.stats.files_failed += 1
            
            return result
            
        except Exception as e:
            result.add_error(f"Directory processing failed: {str(e)}")
            return result
    
    def _process_single_file(self, 
                           csv_file: Path, 
                           output_dir: Path,
                           base_dir: Optional[Path] = None) -> ProcessingResult:
        """Process single CSV file."""
        result = ProcessingResult(success=True)
        
        try:
            self.logger.debug(f"Processing CSV file", file=str(csv_file))
            
            # Determine output filename
            if base_dir:
                # Use parent directory name as output filename
                output_filename = csv_file.parent.name + ".csv"
            else:
                # Use original filename
                output_filename = csv_file.name
            
            output_path = output_dir / output_filename
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Extract and process data
            extracted_data = self._extract_columns_from_file(csv_file)
            
            if not extracted_data:
                result.add_warning(f"No valid data extracted from {csv_file}")
                return result
            
            # Create DataFrame and save
            df = pd.DataFrame(extracted_data, columns=self.column_mapping.target_columns)
            df.to_csv(output_path, index=False)
            
            self.stats.files_processed += 1
            
            result.metadata.update({
                "input_file": str(csv_file),
                "output_file": str(output_path),
                "rows_extracted": len(extracted_data)
            })
            
            self.logger.debug(f"CSV processing completed", 
                            input=str(csv_file),
                            output=str(output_path),
                            rows=len(extracted_data))
            
        except Exception as e:
            result.add_error(f"File processing failed: {str(e)}")
            self.stats.files_failed += 1
        
        return result
    
    def _extract_columns_from_file(self, csv_file: Path) -> List[List[str]]:
        """Extract specified columns from CSV file."""
        extracted_data = []
        
        try:
            with open(csv_file, mode='r', encoding='utf-8', errors='replace') as file:
                csv_reader = csv.reader(file)
                
                # Skip header if configured
                if self.column_mapping.header_row:
                    try:
                        next(csv_reader)
                    except StopIteration:
                        return extracted_data
                
                for row_num, row in enumerate(csv_reader, start=1):
                    self.stats.rows_processed += 1
                    
                    try:
                        # Check minimum column count
                        if len(row) < self.column_mapping.required_columns:
                            self.logger.debug(f"Insufficient columns in row {row_num}", 
                                            expected=self.column_mapping.required_columns,
                                            got=len(row))
                            self.stats.rows_failed += 1
                            continue
                        
                        # Extract specified columns
                        extracted_row = []
                        for col_idx in self.column_mapping.column_indices:
                            if col_idx < len(row):
                                cell_value = row[col_idx]
                                
                                # Apply data cleaning
                                if self.cleaning_config.remove_parentheses or self.cleaning_config.normalize_whitespace:
                                    cell_value = self.data_cleaner.clean_text(cell_value)
                                    self.stats.data_cleaning_applied += 1
                                
                                extracted_row.append(cell_value)
                            else:
                                extracted_row.append("")
                        
                        # Validate extracted data
                        if self.cleaning_config.validate_data_types:
                            is_valid, validation_errors = self.data_cleaner.validate_row_data(
                                extracted_row, len(self.column_mapping.target_columns)
                            )
                            
                            if not is_valid:
                                for error in validation_errors:
                                    self.logger.debug(f"Row validation failed", 
                                                    row=row_num,
                                                    error=error)
                                self.stats.rows_failed += 1
                                continue
                        
                        # Skip empty rows if configured
                        if (self.cleaning_config.remove_empty_rows and 
                            all(not cell.strip() for cell in extracted_row)):
                            continue
                        
                        extracted_data.append(extracted_row)
                        self.stats.rows_extracted += 1
                        
                    except Exception as e:
                        self.logger.warning(f"Failed to process row {row_num}", 
                                          file=str(csv_file),
                                          error=str(e))
                        self.stats.rows_failed += 1
                        continue
            
            return extracted_data
            
        except Exception as e:
            self.logger.error(f"Failed to extract columns from file", 
                            file=str(csv_file),
                            error=str(e))
            raise
    
    def get_processing_statistics(self) -> Dict[str, Any]:
        """Get comprehensive processing statistics."""
        base_stats = self.get_statistics()
        
        success_rate = 0.0
        if self.stats.rows_processed > 0:
            success_rate = self.stats.rows_extracted / self.stats.rows_processed
        
        return {
            **base_stats,
            "files_processed": self.stats.files_processed,
            "files_failed": self.stats.files_failed,
            "rows_processed": self.stats.rows_processed,
            "rows_extracted": self.stats.rows_extracted,
            "rows_failed": self.stats.rows_failed,
            "data_cleaning_applied": self.stats.data_cleaning_applied,
            "extraction_success_rate": success_rate
        }
    
    def validate_input(self, input_data: Path) -> bool:
        """Validate input data."""
        try:
            input_path = Path(input_data) if isinstance(input_data, str) else input_data
            return input_path.exists()
        except Exception:
            return False


def main():
    """Main entry point for column extraction."""
    # Configuration
    input_dir = Path("/groups/itay_mayrose/alongonda/datasets/full_genomes/ensembl/organisms")
    output_dir = Path("/groups/itay_mayrose/alongonda/datasets/full_genomes/ensembl/processed_annotations_test")
    
    # Setup configuration
    config = BioinformaticsConfig(
        input_dir=input_dir,
        output_dir=output_dir,
        max_workers=4,
        log_level="INFO"
    )
    
    column_mapping = ColumnMapping(
        source_columns=["id", "", "", "start", "end", "sequence"],
        target_columns=["id", "start", "end", "sequence"],
        column_indices=[0, 3, 4, 5],
        required_columns=6,
        header_row=True
    )
    
    cleaning_config = DataCleaningConfig(
        remove_parentheses=True,
        normalize_whitespace=True,
        remove_empty_rows=True,
        validate_data_types=True
    )
    
    # Initialize extractor
    try:
        extractor = ColumnExtractor(
            config=config,
            column_mapping=column_mapping,
            cleaning_config=cleaning_config
        )
        
        # Run extraction
        result = extractor.run(input_dir, output_dir=output_dir)
        
        # Display results
        if result.success:
            print("✅ Column extraction completed successfully!")
            stats = extractor.get_processing_statistics()
            print(f"   Files processed: {stats['files_processed']}")
            print(f"   Rows extracted: {stats['rows_extracted']}")
            print(f"   Rows failed: {stats['rows_failed']}")
            print(f"   Success rate: {stats['extraction_success_rate']:.2%}")
            print(f"   Processing time: {result.processing_time:.2f} seconds")
        else:
            print("❌ Column extraction failed!")
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