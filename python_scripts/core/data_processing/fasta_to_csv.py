"""
Industry-Level FASTA to CSV Converter.

This module provides object-oriented FASTA to CSV conversion with comprehensive
header parsing, parallel processing, and robust error handling.
"""

import os
from typing import Dict, List, Optional, Any, Union, Tuple
from dataclasses import dataclass, field
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import re

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..base_classes import (
    BaseProcessor,
    BioinformaticsConfig,
    ProcessingResult,
    SequenceProcessor,
    ParallelProcessor
)


@dataclass
class FASTAConversionConfig:
    """Configuration for FASTA to CSV conversion."""
    
    header_formats: List[str] = field(default_factory=lambda: [
        "mgc_candidate|gene_id|kegg_id|source_file",
        "gene_id|mgc_id|start|end",
        "identifier|annotation|start|end",
        "locus_tag|organism|start|end"
    ])
    delimiter: str = "|"
    sequence_column_name: str = "sequence"
    auto_detect_format: bool = True
    validate_sequences: bool = True
    min_sequence_length: int = 1
    include_sequence_stats: bool = False
    
    def get_format_patterns(self) -> List[re.Pattern]:
        """Get compiled regex patterns for header formats."""
        patterns = []
        for fmt in self.header_formats:
            # Convert format to regex pattern
            pattern = fmt.replace("|", r"\|").replace("start", r"\d+").replace("end", r"\d+")
            patterns.append(re.compile(pattern))
        return patterns


@dataclass
class FASTARecord:
    """Parsed FASTA record with structured data."""
    
    raw_header: str
    sequence: str
    parsed_fields: Dict[str, Any] = field(default_factory=dict)
    format_type: Optional[str] = None
    sequence_length: int = 0
    gc_content: float = 0.0
    
    def __post_init__(self):
        """Calculate sequence statistics."""
        self.sequence_length = len(self.sequence)
        if self.sequence:
            gc_count = self.sequence.upper().count('G') + self.sequence.upper().count('C')
            self.gc_content = gc_count / len(self.sequence) if len(self.sequence) > 0 else 0.0
    
    def to_dict(self, include_stats: bool = False) -> Dict[str, Any]:
        """Convert to dictionary for DataFrame creation."""
        result = dict(self.parsed_fields)
        result['sequence'] = self.sequence
        
        if include_stats:
            result['sequence_length'] = self.sequence_length
            result['gc_content'] = self.gc_content
        
        return result


class FASTAHeaderParser:
    """Parser for different FASTA header formats."""
    
    def __init__(self, config: FASTAConversionConfig, logger):
        """Initialize header parser."""
        self.config = config
        self.logger = logger
        self.format_stats = {fmt: 0 for fmt in config.header_formats}
    
    def parse_header(self, header: str) -> Tuple[Dict[str, Any], Optional[str]]:
        """
        Parse FASTA header into structured fields.
        
        Args:
            header: FASTA header line (without >)
            
        Returns:
            Tuple of (parsed_fields, format_type)
        """
        try:
            # Clean header
            header = header.strip()
            
            # Split by delimiter
            parts = [part.strip() for part in header.split(self.config.delimiter)]
            
            if len(parts) < 2:
                # Single field header - use as identifier
                return {"identifier": header}, "single_field"
            
            # Try to match known formats
            if self.config.auto_detect_format:
                parsed_fields, format_type = self._auto_detect_format(parts)
                if parsed_fields:
                    self.format_stats[format_type] += 1
                    return parsed_fields, format_type
            
            # Try each known format
            for format_str in self.config.header_formats:
                format_fields = format_str.split(self.config.delimiter)
                
                if len(parts) == len(format_fields):
                    parsed_fields = self._parse_with_format(parts, format_fields)
                    if parsed_fields:
                        self.format_stats[format_str] += 1
                        return parsed_fields, format_str
            
            # Fallback - create generic fields
            generic_fields = {f"field_{i}": part for i, part in enumerate(parts)}
            return generic_fields, "generic"
            
        except Exception as e:
            self.logger.warning(f"Failed to parse header", header=header, error=str(e))
            return {"identifier": header}, "unparseable"
    
    def _auto_detect_format(self, parts: List[str]) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
        """Auto-detect header format based on content."""
        try:
            if len(parts) == 4:
                # Check if last two parts are integers (coordinates)
                try:
                    start = int(parts[2])
                    end = int(parts[3])
                    
                    # Coordinate format: gene_id | mgc_id | start | end
                    return {
                        'gene_id': parts[0],
                        'mgc_id': parts[1],
                        'start': start,
                        'end': end
                    }, "gene_id|mgc_id|start|end"
                    
                except ValueError:
                    # String format: mgc_candidate | gene_id | kegg_id | source_file
                    return {
                        'mgc_candidate': parts[0],
                        'gene_id': parts[1],
                        'kegg_id': parts[2],
                        'source_file': parts[3]
                    }, "mgc_candidate|gene_id|kegg_id|source_file"
            
            elif len(parts) == 3:
                # Three field format
                return {
                    'identifier': parts[0],
                    'annotation': parts[1],
                    'source': parts[2]
                }, "identifier|annotation|source"
            
            elif len(parts) == 2:
                # Two field format
                return {
                    'identifier': parts[0],
                    'description': parts[1]
                }, "identifier|description"
            
            return None, None
            
        except Exception:
            return None, None
    
    def _parse_with_format(self, parts: List[str], format_fields: List[str]) -> Optional[Dict[str, Any]]:
        """Parse parts using specific format."""
        try:
            parsed = {}
            
            for i, field_name in enumerate(format_fields):
                if i < len(parts):
                    value = parts[i]
                    
                    # Try to convert coordinates to integers
                    if field_name in ['start', 'end'] and value.isdigit():
                        parsed[field_name] = int(value)
                    else:
                        parsed[field_name] = value
            
            return parsed
            
        except Exception:
            return None
    
    def get_format_statistics(self) -> Dict[str, int]:
        """Get statistics on header formats encountered."""
        return dict(self.format_stats)


class FASTASequenceValidator:
    """Validator for FASTA sequences."""
    
    def __init__(self, config: FASTAConversionConfig, logger):
        """Initialize sequence validator."""
        self.config = config
        self.logger = logger
    
    def validate_sequence(self, sequence: str, record_id: str) -> Tuple[bool, List[str]]:
        """
        Validate FASTA sequence.
        
        Args:
            sequence: Sequence to validate
            record_id: Record identifier for error reporting
            
        Returns:
            Tuple of (is_valid, error_messages)
        """
        errors = []
        
        if not sequence:
            errors.append(f"Empty sequence for {record_id}")
            return False, errors
        
        # Check minimum length
        if len(sequence) < self.config.min_sequence_length:
            errors.append(f"Sequence too short for {record_id}: {len(sequence)} < {self.config.min_sequence_length}")
        
        # Check for valid characters (protein or nucleotide)
        sequence_upper = sequence.upper()
        
        # Check if nucleotide
        nucleotide_chars = set("ATCGNU-")
        if set(sequence_upper).issubset(nucleotide_chars):
            # Valid nucleotide sequence
            pass
        else:
            # Check if protein
            protein_chars = set("ACDEFGHIKLMNPQRSTVWY*-")
            invalid_chars = set(sequence_upper) - protein_chars
            if invalid_chars:
                errors.append(f"Invalid characters in sequence for {record_id}: {invalid_chars}")
        
        return len(errors) == 0, errors


class FASTAToCSVConverter(SequenceProcessor):
    """
    Professional FASTA to CSV converter with comprehensive functionality.
    """
    
    def __init__(self, 
                 config: BioinformaticsConfig,
                 conversion_config: FASTAConversionConfig):
        """
        Initialize FASTA to CSV converter.
        
        Args:
            config: Base processing configuration
            conversion_config: Conversion-specific configuration
        """
        super().__init__(config, name="FASTAToCSVConverter")
        
        self.conversion_config = conversion_config
        self.header_parser = FASTAHeaderParser(conversion_config, self.logger)
        self.sequence_validator = FASTASequenceValidator(conversion_config, self.logger)
        
        # Statistics
        self.conversion_stats = {
            "files_processed": 0,
            "files_failed": 0,
            "records_converted": 0,
            "records_failed": 0,
            "validation_errors": 0,
        }
    
    def process(self, input_data: Union[Path, str], **kwargs) -> ProcessingResult:
        """
        Process FASTA to CSV conversion.
        
        Args:
            input_data: Path to FASTA file, directory, or file list
            **kwargs: Additional processing parameters
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            input_path = Path(input_data) if isinstance(input_data, str) else input_data
            
            self.logger.info(f"Starting FASTA to CSV conversion", input=str(input_path))
            
            # Determine input type and process accordingly
            if input_path.is_file():
                if input_path.suffix.lower() in ['.fa', '.fasta', '.fas']:
                    # Single FASTA file
                    conversion_result = self._convert_single_file(input_path)
                    result = result.merge(conversion_result)
                elif input_path.suffix.lower() == '.txt':
                    # File list
                    list_result = self._process_file_list(input_path)
                    result = result.merge(list_result)
                else:
                    result.add_error(f"Unsupported file type: {input_path}")
            
            elif input_path.is_dir():
                # Directory of FASTA files
                directory_result = self._convert_directory(input_path)
                result = result.merge(directory_result)
            
            else:
                result.add_error(f"Input path does not exist: {input_path}")
            
            # Add format statistics to results
            format_stats = self.header_parser.get_format_statistics()
            
            result.metadata.update({
                "conversion_stats": self.conversion_stats,
                "format_statistics": format_stats,
                "input_path": str(input_path),
                "output_directory": str(self.config.output_dir)
            })
            
            self.logger.info(f"FASTA to CSV conversion completed",
                           files_processed=self.conversion_stats["files_processed"],
                           records_converted=self.conversion_stats["records_converted"])
            
        except Exception as e:
            result.add_error(f"Conversion failed: {str(e)}")
        
        return result
    
    def _process_file_list(self, list_file: Path) -> ProcessingResult:
        """
        Process FASTA files from a text file list.
        
        Args:
            list_file: Path to text file containing FASTA file paths
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            # Read file list
            with open(list_file, 'r') as f:
                fasta_files = [line.strip() for line in f if line.strip()]
            
            if not fasta_files:
                result.add_warning("No FASTA files found in list")
                return result
            
            self.logger.info(f"Found FASTA files in list", count=len(fasta_files))
            
            # Convert paths to Path objects and validate
            valid_files = []
            for fasta_file in fasta_files:
                fasta_path = Path(fasta_file)
                if fasta_path.exists():
                    valid_files.append(fasta_path)
                else:
                    result.add_warning(f"FASTA file not found: {fasta_path}")
            
            # Process files in parallel
            if self.config.max_workers == 1:
                # Sequential processing
                for fasta_file in valid_files:
                    file_result = self._convert_single_file(fasta_file)
                    if not file_result.success:
                        result.errors.extend(file_result.errors)
                        result.warnings.extend(file_result.warnings)
            else:
                # Parallel processing
                with ProcessPoolExecutor(max_workers=self.config.max_workers) as executor:
                    future_to_file = {
                        executor.submit(self._convert_single_file, fasta_file): fasta_file
                        for fasta_file in valid_files
                    }
                    
                    for future in as_completed(future_to_file):
                        fasta_file = future_to_file[future]
                        try:
                            file_result = future.result(timeout=self.config.timeout)
                            if not file_result.success:
                                result.errors.extend(file_result.errors)
                                result.warnings.extend(file_result.warnings)
                        except Exception as e:
                            result.add_error(f"Failed to process {fasta_file}: {str(e)}")
                            self.conversion_stats["files_failed"] += 1
            
            return result
            
        except Exception as e:
            result.add_error(f"File list processing failed: {str(e)}")
            return result
    
    def _convert_directory(self, input_dir: Path) -> ProcessingResult:
        """Convert all FASTA files in directory."""
        result = ProcessingResult(success=True)
        
        try:
            # Find FASTA files
            fasta_extensions = ['.fa', '.fasta', '.fas']
            fasta_files = []
            
            for ext in fasta_extensions:
                fasta_files.extend(input_dir.rglob(f"*{ext}"))
            
            if not fasta_files:
                result.add_warning("No FASTA files found in directory")
                return result
            
            self.logger.info(f"Found FASTA files", count=len(fasta_files))
            
            # Process files
            for fasta_file in fasta_files:
                file_result = self._convert_single_file(fasta_file)
                if not file_result.success:
                    result.errors.extend(file_result.errors)
                    result.warnings.extend(file_result.warnings)
            
            return result
            
        except Exception as e:
            result.add_error(f"Directory conversion failed: {str(e)}")
            return result
    
    def _convert_single_file(self, fasta_file: Path) -> ProcessingResult:
        """
        Convert single FASTA file to CSV.
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            self.logger.debug(f"Converting FASTA file", file=str(fasta_file))
            
            # Parse FASTA file
            records = self._parse_fasta_file(fasta_file)
            
            if not records:
                result.add_warning(f"No valid records found in {fasta_file}")
                return result
            
            # Convert to DataFrame
            df = self._records_to_dataframe(records)
            
            # Determine output path
            output_path = self.config.output_dir / f"{fasta_file.stem}.csv"
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Write CSV
            df.to_csv(output_path, index=False)
            
            # Update statistics
            self.conversion_stats["files_processed"] += 1
            self.conversion_stats["records_converted"] += len(records)
            
            result.metadata.update({
                "input_file": str(fasta_file),
                "output_file": str(output_path),
                "records_converted": len(records)
            })
            
            self.logger.debug(f"FASTA conversion completed", 
                            input=str(fasta_file),
                            output=str(output_path),
                            records=len(records))
            
        except Exception as e:
            result.add_error(f"File conversion failed: {str(e)}")
            self.conversion_stats["files_failed"] += 1
        
        return result
    
    def _parse_fasta_file(self, fasta_file: Path) -> List[FASTARecord]:
        """
        Parse FASTA file into structured records.
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            List of FASTARecord objects
        """
        records = []
        
        try:
            for seq_record in SeqIO.parse(fasta_file, "fasta"):
                # Parse header
                parsed_fields, format_type = self.header_parser.parse_header(seq_record.description)
                
                # Get sequence
                sequence = str(seq_record.seq)
                
                # Validate sequence if configured
                if self.conversion_config.validate_sequences:
                    is_valid, validation_errors = self.sequence_validator.validate_sequence(
                        sequence, seq_record.id
                    )
                    
                    if not is_valid:
                        for error in validation_errors:
                            self.logger.warning(error, file=str(fasta_file))
                        self.conversion_stats["validation_errors"] += len(validation_errors)
                        self.conversion_stats["records_failed"] += 1
                        continue
                
                # Create record
                record = FASTARecord(
                    raw_header=seq_record.description,
                    sequence=sequence,
                    parsed_fields=parsed_fields,
                    format_type=format_type
                )
                
                records.append(record)
            
            return records
            
        except Exception as e:
            self.logger.error(f"Failed to parse FASTA file", 
                            file=str(fasta_file),
                            error=str(e))
            raise
    
    def _records_to_dataframe(self, records: List[FASTARecord]) -> pd.DataFrame:
        """Convert FASTA records to DataFrame."""
        try:
            data = []
            
            for record in records:
                record_dict = record.to_dict(self.conversion_config.include_sequence_stats)
                data.append(record_dict)
            
            return pd.DataFrame(data)
            
        except Exception as e:
            self.logger.error(f"Failed to convert records to DataFrame", error=str(e))
            raise
    
    def get_conversion_statistics(self) -> Dict[str, Any]:
        """Get comprehensive conversion statistics."""
        base_stats = self.get_statistics()
        format_stats = self.header_parser.get_format_statistics()
        
        success_rate = 0.0
        if self.conversion_stats["files_processed"] + self.conversion_stats["files_failed"] > 0:
            success_rate = (self.conversion_stats["files_processed"] / 
                          (self.conversion_stats["files_processed"] + self.conversion_stats["files_failed"]))
        
        return {
            **base_stats,
            **self.conversion_stats,
            "format_statistics": format_stats,
            "success_rate": success_rate,
            "records_per_file": (
                self.conversion_stats["records_converted"] / 
                max(1, self.conversion_stats["files_processed"])
            )
        }
    
    def validate_input(self, input_data: Union[Path, str]) -> bool:
        """Validate input data."""
        try:
            input_path = Path(input_data) if isinstance(input_data, str) else input_data
            return input_path.exists()
        except Exception:
            return False


def main():
    """Main entry point for FASTA to CSV conversion."""
    # Configuration
    fasta_list_path = Path('/groups/itay_mayrose/alongonda/Plant_MGC/kegg_metabolic_output_g3_slurm_no_chloroplast/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/merged_list.txt')
    output_folder = Path('/groups/itay_mayrose/alongonda/Plant_MGC/kegg_metabolic_output_g3_slurm_no_chloroplast/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mgc_candidates_csv_files')
    
    # Setup configuration
    config = BioinformaticsConfig(
        input_dir=fasta_list_path.parent,
        output_dir=output_folder,
        max_workers=30,
        log_level="INFO"
    )
    
    conversion_config = FASTAConversionConfig(
        auto_detect_format=True,
        validate_sequences=True,
        include_sequence_stats=True,
        min_sequence_length=10
    )
    
    # Initialize converter
    try:
        converter = FASTAToCSVConverter(
            config=config,
            conversion_config=conversion_config
        )
        
        # Run conversion
        result = converter.run(fasta_list_path)
        
        # Display results
        if result.success:
            print("✅ FASTA to CSV conversion completed successfully!")
            stats = converter.get_conversion_statistics()
            print(f"   Files processed: {stats['files_processed']}")
            print(f"   Records converted: {stats['records_converted']}")
            print(f"   Success rate: {stats['success_rate']:.2%}")
            print(f"   Processing time: {result.processing_time:.2f} seconds")
            
            # Display format statistics
            format_stats = stats['format_statistics']
            print(f"\nHeader Format Statistics:")
            for fmt, count in format_stats.items():
                if count > 0:
                    print(f"   {fmt}: {count} records")
        else:
            print("❌ FASTA to CSV conversion failed!")
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
