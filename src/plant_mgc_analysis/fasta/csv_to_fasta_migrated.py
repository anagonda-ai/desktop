"""
Industry-Level CSV to FASTA Converter.

This module provides object-oriented CSV to FASTA conversion with comprehensive
error handling, parallel processing, and performance optimization.
"""

from ..config.settings import get_settings
from ..core.base import BioinformaticsProcessor
from ..core.exceptions import ProcessingError, ValidationError
from ..utils.file_operations import file_manager
from loguru import logger
import csv
import sys
import os
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from ..base_classes import (
    BaseProcessor,
    BioinformaticsConfig,
    ProcessingResult,
    SequenceProcessor,
    FileSystemProcessor
)


@dataclass
class CSVConversionConfig:
    """Configuration for CSV to FASTA conversion."""
    
    sequence_column: str = "Translation"
    locus_tag_column: str = "Locus_Tag"
    protein_id_column: str = "Protein_ID"
    gene_column: str = "Gene"
    start_column: str = "Start"
    end_column: str = "End"
    header_format: str = "{id} | {source} | {start} | {end}"
    preserve_directory_structure: bool = True
    validate_sequences: bool = True
    min_sequence_length: int = 1
    
    def get_identifier_columns(self) -> List[str]:
        """Get list of identifier columns in priority order."""
        return [self.locus_tag_column, self.protein_id_column, self.gene_column]


@dataclass
class ConversionRecord:
    """Record for a single CSV to FASTA conversion."""
    
    identifier: str
    sequence: str
    source_file: str
    start_position: Optional[str] = None
    end_position: Optional[str] = None
    additional_info: Dict[str, Any] = None
    
    def to_seqrecord(self, header_format: str) -> SeqRecord:
        """Convert to BioPython SeqRecord."""
        header = header_format.format(
            id=self.identifier,
            source=self.source_file,
            start=self.start_position or "unknown",
            end=self.end_position or "unknown"
        )
        
        return SeqRecord(
            Seq(self.sequence),
            id=self.identifier,
            description=header
        )


class CSVValidator:
    """Validator for CSV data integrity."""
    
    def __init__(self, config: CSVConversionConfig, logger):
        """Initialize CSV validator."""
        self.config = config
        self.logger = logger
    
    def validate_csv_structure(self, csv_path: Path) -> Tuple[bool, List[str]]:
        """
        Validate CSV file structure and required columns.
        
        Args:
            csv_path: Path to CSV file
            
        Returns:
            Tuple of (is_valid, error_messages)
        """
        errors = []
        
        try:
            # Read first few rows to check structure
            df = pd.read_csv(csv_path, nrows=5)
            
            # Check for required sequence column
            if self.config.sequence_column not in df.columns:
                errors.append(f"Missing required sequence column: {self.config.sequence_column}")
            
            # Check for at least one identifier column
            id_columns = self.config.get_identifier_columns()
            if not any(col in df.columns for col in id_columns):
                errors.append(f"Missing identifier columns. Need at least one of: {id_columns}")
            
            # Check for position columns if configured
            position_columns = [self.config.start_column, self.config.end_column]
            missing_positions = [col for col in position_columns if col not in df.columns]
            if missing_positions:
                self.logger.warning(f"Missing position columns", 
                                  file=str(csv_path),
                                  missing=missing_positions)
            
            return len(errors) == 0, errors
            
        except Exception as e:
            errors.append(f"Failed to read CSV file: {str(e)}")
            return False, errors
    
    def validate_sequence(self, sequence: str, row_index: int) -> Tuple[bool, List[str]]:
        """
        Validate individual sequence.
        
        Args:
            sequence: Protein sequence to validate
            row_index: Row number for error reporting
            
        Returns:
            Tuple of (is_valid, error_messages)
        """
        errors = []
        
        if not sequence or pd.isna(sequence):
            errors.append(f"Row {row_index}: Empty sequence")
            return False, errors
        
        sequence_str = str(sequence).strip()
        
        # Check minimum length
        if len(sequence_str) < self.config.min_sequence_length:
            errors.append(f"Row {row_index}: Sequence too short ({len(sequence_str)} < {self.config.min_sequence_length})")
        
        # Check for valid protein characters
        valid_aa_chars = set("ACDEFGHIKLMNPQRSTVWY*-")
        invalid_chars = set(sequence_str.upper()) - valid_aa_chars
        if invalid_chars:
            errors.append(f"Row {row_index}: Invalid amino acid characters: {invalid_chars}")
        
        return len(errors) == 0, errors


class CSVToFASTAConverter(SequenceProcessor):
    """
    Professional CSV to FASTA converter with comprehensive functionality.
    """
    
    def __init__(self, 
                 config: BioinformaticsConfig,
                 conversion_config: CSVConversionConfig):
        """
        Initialize CSV to FASTA converter.
        
        Args:
            config: Base processing configuration
            conversion_config: Conversion-specific configuration
        """
        super().__init__(config, name="CSVToFASTAConverter")
        
        self.conversion_config = conversion_config
        self.validator = CSVValidator(conversion_config, self.logger)
        
        # Statistics
        self.conversion_stats = {
            "files_processed": 0,
            "files_failed": 0,
            "records_converted": 0,
            "records_failed": 0,
            "validation_errors": 0,
        }
    
    def process(self, input_data: Path, **kwargs) -> ProcessingResult:
        """
        Process CSV to FASTA conversion.
        
        Args:
            input_data: Path to input directory or single CSV file
            **kwargs: Additional processing parameters
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            self.logger.info(f"Starting CSV to FASTA conversion", input=str(input_data))
            
            # Determine if input is file or directory
            if input_data.is_file():
                if input_data.suffix.lower() == '.csv':
                    conversion_result = self._convert_single_file(input_data)
                    result = result.merge(conversion_result)
                else:
                    result.add_error(f"Input file is not a CSV: {input_data}")
            
            elif input_data.is_dir():
                directory_result = self._convert_directory(input_data)
                result = result.merge(directory_result)
            
            else:
                result.add_error(f"Input path does not exist: {input_data}")
            
            result.metadata.update({
                "conversion_stats": self.conversion_stats,
                "input_path": str(input_data),
                "output_directory": str(self.config.output_dir)
            })
            
            self.logger.info(f"CSV to FASTA conversion completed",
                           files_processed=self.conversion_stats["files_processed"],
                           records_converted=self.conversion_stats["records_converted"])
            
        except Exception as e:
            result.add_error(f"Conversion failed: {str(e)}")
        
        return result
    
    def _convert_directory(self, input_dir: Path) -> ProcessingResult:
        """
        Convert all CSV files in directory.
        
        Args:
            input_dir: Input directory path
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            # Find all CSV files
            csv_files = list(input_dir.rglob("*.csv"))
            
            if not csv_files:
                result.add_warning("No CSV files found in directory")
                return result
            
            self.logger.info(f"Found CSV files", count=len(csv_files))
            
            # Process files in parallel
            if self.config.max_workers == 1:
                # Sequential processing
                for csv_file in csv_files:
                    file_result = self._convert_single_file(csv_file, input_dir)
                    if not file_result.success:
                        result.errors.extend(file_result.errors)
                        result.warnings.extend(file_result.warnings)
            else:
                # Parallel processing
                with ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
                    future_to_file = {
                        executor.submit(self._convert_single_file, csv_file, input_dir): csv_file
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
                            self.conversion_stats["files_failed"] += 1
            
            return result
            
        except Exception as e:
            result.add_error(f"Directory conversion failed: {str(e)}")
            return result
    
    def _convert_single_file(self, 
                           csv_file: Path, 
                           base_dir: Optional[Path] = None) -> ProcessingResult:
        """
        Convert single CSV file to FASTA.
        
        Args:
            csv_file: Path to CSV file
            base_dir: Base directory for relative path calculation
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            self.logger.debug(f"Converting CSV file", file=str(csv_file))
            
            # Validate CSV structure
            is_valid, validation_errors = self.validator.validate_csv_structure(csv_file)
            if not is_valid:
                for error in validation_errors:
                    result.add_error(error)
                self.conversion_stats["files_failed"] += 1
                return result
            
            # Determine output path
            output_path = self._get_output_path(csv_file, base_dir)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Process CSV data
            records = self._process_csv_data(csv_file)
            
            if not records:
                result.add_warning(f"No valid records found in {csv_file}")
                return result
            
            # Write FASTA file
            self._write_fasta_file(records, output_path)
            
            # Update statistics
            self.conversion_stats["files_processed"] += 1
            self.conversion_stats["records_converted"] += len(records)
            
            result.metadata.update({
                "input_file": str(csv_file),
                "output_file": str(output_path),
                "records_converted": len(records)
            })
            
            self.logger.debug(f"CSV conversion completed", 
                            input=str(csv_file),
                            output=str(output_path),
                            records=len(records))
            
        except Exception as e:
            result.add_error(f"File conversion failed: {str(e)}")
            self.conversion_stats["files_failed"] += 1
        
        return result
    
    def _process_csv_data(self, csv_file: Path) -> List[ConversionRecord]:
        """
        Process CSV data and extract conversion records.
        
        Args:
            csv_file: Path to CSV file
            
        Returns:
            List of ConversionRecord objects
        """
        records = []
        
        try:
            df = pd.read_csv(csv_file)
            source_name = csv_file.stem
            
            for index, row in df.iterrows():
                try:
                    # Extract identifier (first available from priority list)
                    identifier = None
                    for col in self.conversion_config.get_identifier_columns():
                        if col in df.columns and pd.notna(row[col]):
                            identifier = str(row[col])
                            break
                    
                    if not identifier:
                        self.logger.warning(f"No identifier found for row {index}", file=str(csv_file))
                        self.conversion_stats["records_failed"] += 1
                        continue
                    
                    # Extract sequence
                    sequence = row.get(self.conversion_config.sequence_column, "")
                    if pd.isna(sequence):
                        sequence = ""
                    
                    # Validate sequence if configured
                    if self.conversion_config.validate_sequences:
                        is_valid, validation_errors = self.validator.validate_sequence(sequence, index)
                        if not is_valid:
                            for error in validation_errors:
                                self.logger.warning(error, file=str(csv_file))
                            self.conversion_stats["validation_errors"] += len(validation_errors)
                            self.conversion_stats["records_failed"] += 1
                            continue
                    
                    # Extract position information
                    start_pos = row.get(self.conversion_config.start_column)
                    end_pos = row.get(self.conversion_config.end_column)
                    
                    # Create conversion record
                    record = ConversionRecord(
                        identifier=identifier,
                        sequence=str(sequence).strip(),
                        source_file=source_name,
                        start_position=str(start_pos) if pd.notna(start_pos) else None,
                        end_position=str(end_pos) if pd.notna(end_pos) else None,
                        additional_info=row.to_dict()
                    )
                    
                    records.append(record)
                    
                except Exception as e:
                    self.logger.warning(f"Failed to process row {index}", 
                                      file=str(csv_file),
                                      error=str(e))
                    self.conversion_stats["records_failed"] += 1
            
            return records
            
        except Exception as e:
            self.logger.error(f"Failed to process CSV data", 
                            file=str(csv_file),
                            error=str(e))
            raise
    
    def _write_fasta_file(self, records: List[ConversionRecord], output_path: Path) -> None:
        """
        Write conversion records to FASTA file.
        
        Args:
            records: List of conversion records
            output_path: Output FASTA file path
        """
        try:
            seq_records = []
            
            for record in records:
                seq_record = record.to_seqrecord(self.conversion_config.header_format)
                seq_records.append(seq_record)
            
            with open(output_path, 'w') as f:
                SeqIO.write(seq_records, f, "fasta")
            
        except Exception as e:
            self.logger.error(f"Failed to write FASTA file", 
                            file=str(output_path),
                            error=str(e))
            raise
    
    def _get_output_path(self, csv_file: Path, base_dir: Optional[Path] = None) -> Path:
        """
        Determine output path for FASTA file.
        
        Args:
            csv_file: Input CSV file path
            base_dir: Base directory for relative path calculation
            
        Returns:
            Output FASTA file path
        """
        if self.conversion_config.preserve_directory_structure and base_dir:
            # Preserve directory structure
            relative_path = csv_file.relative_to(base_dir)
            output_dir = self.config.output_dir / relative_path.parent
        else:
            # Flat output structure
            output_dir = self.config.output_dir
        
        fasta_filename = f"{csv_file.stem}.fasta"
        return output_dir / fasta_filename
    
    def get_conversion_statistics(self) -> Dict[str, Any]:
        """Get comprehensive conversion statistics."""
        base_stats = self.get_statistics()
        
        success_rate = 0.0
        if self.conversion_stats["files_processed"] + self.conversion_stats["files_failed"] > 0:
            success_rate = (self.conversion_stats["files_processed"] / 
                          (self.conversion_stats["files_processed"] + self.conversion_stats["files_failed"]))
        
        return {
            **base_stats,
            **self.conversion_stats,
            "success_rate": success_rate,
            "records_per_file": (
                self.conversion_stats["records_converted"] / 
                max(1, self.conversion_stats["files_processed"])
            )
        }
    
    def validate_input(self, input_data: Path) -> bool:
        """
        Validate input data.
        
        Args:
            input_data: Input path to validate
            
        Returns:
            True if valid
        """
        if not isinstance(input_data, Path):
            input_data = Path(input_data)
        
        return input_data.exists()


def main():
    """Main entry point for CSV to FASTA conversion."""
    # Configuration
    input_dir = Path("/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/csv_files")
    output_dir = Path("/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/fasta_files")
    
    # Setup configuration
    config = BioinformaticsConfig(
        input_dir=input_dir,
        output_dir=output_dir,
        max_workers=8,
        log_level="INFO"
    )
    
    conversion_config = CSVConversionConfig(
        sequence_column="Translation",
        preserve_directory_structure=True,
        validate_sequences=True,
        min_sequence_length=10
    )
    
    # Initialize converter
    try:
        converter = CSVToFASTAConverter(
            config=config,
            conversion_config=conversion_config
        )
        
        # Run conversion
        result = converter.run(input_dir)
        
        # Display results
        if result.success:
            logger.info("✅ CSV to FASTA conversion completed successfully!")
            stats = converter.get_conversion_statistics()
            logger.info(f"   Files processed: {stats['files_processed']}")
            logger.info(f"   Records converted: {stats['records_converted']}")
            logger.info(f"   Success rate: {stats['success_rate']:.2%}")
            logger.info(f"   Processing time: {result.processing_time:.2f} seconds")
        else:
            logger.info("❌ CSV to FASTA conversion failed!")
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