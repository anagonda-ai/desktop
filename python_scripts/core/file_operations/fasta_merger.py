"""
Industry-Level FASTA File Merger.

This module provides object-oriented FASTA file merging with comprehensive
error handling, validation, sequence deduplication, and parallel processing.
"""

import os
from typing import Dict, List, Optional, Any, Set, Tuple
from dataclasses import dataclass, field
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from ..base_classes import (
    BaseProcessor,
    BioinformaticsConfig,
    ProcessingResult,
    SequenceProcessor,
    FileSystemProcessor
)


@dataclass
class FASTAMergeConfig:
    """Configuration for FASTA merging operations."""
    
    add_source_suffix: bool = True
    source_separator: str = "$"
    deduplicate_sequences: bool = False
    deduplicate_by: str = "sequence"  # "sequence", "header", "both"
    validate_sequences: bool = True
    preserve_descriptions: bool = True
    max_file_size_mb: float = 500.0
    chunk_processing: bool = False
    
    def get_source_suffix(self, file_path: Path) -> str:
        """Get source suffix from file path."""
        return file_path.stem


@dataclass
class FASTAFileInfo:
    """Information about a FASTA file to be merged."""
    
    file_path: Path
    source_name: str
    sequence_count: int = 0
    file_size_mb: float = 0.0
    valid_sequences: int = 0
    invalid_sequences: int = 0
    duplicate_sequences: int = 0
    
    def __post_init__(self):
        """Calculate file statistics."""
        if self.file_path.exists():
            self.file_size_mb = self.file_path.stat().st_size / (1024 * 1024)


class SequenceDeduplicator:
    """Handles sequence deduplication logic."""
    
    def __init__(self, config: FASTAMergeConfig, logger):
        """Initialize sequence deduplicator."""
        self.config = config
        self.logger = logger
        self.seen_sequences: Set[str] = set()
        self.seen_headers: Set[str] = set()
        self.duplicate_count = 0
    
    def is_duplicate(self, sequence_record: SeqRecord) -> bool:
        """
        Check if sequence record is a duplicate.
        
        Args:
            sequence_record: BioPython SeqRecord
            
        Returns:
            True if duplicate, False otherwise
        """
        if not self.config.deduplicate_sequences:
            return False
        
        sequence_str = str(sequence_record.seq).upper()
        header_str = sequence_record.id
        
        is_dup = False
        
        if self.config.deduplicate_by == "sequence":
            is_dup = sequence_str in self.seen_sequences
            if not is_dup:
                self.seen_sequences.add(sequence_str)
        elif self.config.deduplicate_by == "header":
            is_dup = header_str in self.seen_headers
            if not is_dup:
                self.seen_headers.add(header_str)
        elif self.config.deduplicate_by == "both":
            is_dup = (sequence_str in self.seen_sequences or 
                     header_str in self.seen_headers)
            if not is_dup:
                self.seen_sequences.add(sequence_str)
                self.seen_headers.add(header_str)
        
        if is_dup:
            self.duplicate_count += 1
            self.logger.debug(f"Duplicate sequence found", 
                            header=header_str,
                            deduplicate_by=self.config.deduplicate_by)
        
        return is_dup
    
    def get_duplicate_count(self) -> int:
        """Get total number of duplicates found."""
        return self.duplicate_count


class FASTASequenceValidator:
    """Validator for FASTA sequences."""
    
    def __init__(self, config: FASTAMergeConfig, logger):
        """Initialize sequence validator."""
        self.config = config
        self.logger = logger
    
    def validate_sequence(self, sequence_record: SeqRecord) -> Tuple[bool, List[str]]:
        """
        Validate FASTA sequence record.
        
        Args:
            sequence_record: BioPython SeqRecord
            
        Returns:
            Tuple of (is_valid, error_messages)
        """
        if not self.config.validate_sequences:
            return True, []
        
        errors = []
        
        # Check for empty sequence
        if len(sequence_record.seq) == 0:
            errors.append(f"Empty sequence for {sequence_record.id}")
        
        # Check for valid characters (basic validation)
        sequence_str = str(sequence_record.seq).upper()
        
        # Check for common invalid characters
        invalid_chars = set(sequence_str) - set("ACDEFGHIKLMNPQRSTVWY*-ATCGUN")
        if invalid_chars:
            # Allow if it's mostly valid (might be mixed nucleotide/protein)
            invalid_ratio = len(invalid_chars) / len(set(sequence_str)) if sequence_str else 0
            if invalid_ratio > 0.1:  # More than 10% invalid characters
                errors.append(f"Invalid characters in sequence {sequence_record.id}: {invalid_chars}")
        
        return len(errors) == 0, errors


class FASTAFileMerger(SequenceProcessor):
    """
    Professional FASTA file merger with comprehensive functionality.
    """
    
    def __init__(self, 
                 config: BioinformaticsConfig,
                 merge_config: FASTAMergeConfig):
        """
        Initialize FASTA file merger.
        
        Args:
            config: Base processing configuration
            merge_config: Merge-specific configuration
        """
        super().__init__(config, name="FASTAFileMerger")
        
        self.merge_config = merge_config
        self.sequence_validator = FASTASequenceValidator(merge_config, self.logger)
        self.deduplicator = SequenceDeduplicator(merge_config, self.logger)
        
        # Statistics
        self.merge_stats = {
            "files_found": 0,
            "files_processed": 0,
            "files_failed": 0,
            "total_sequences_processed": 0,
            "total_sequences_merged": 0,
            "sequences_deduplicated": 0,
            "sequences_failed_validation": 0,
            "total_size_mb": 0.0,
        }
    
    def process(self, input_data: Path, **kwargs) -> ProcessingResult:
        """
        Process FASTA file merging.
        
        Args:
            input_data: Path to directory containing FASTA files
            **kwargs: Additional processing parameters including:
                - output_file: Output FASTA file path
                - output_name: Output file name (without extension)
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            input_dir = Path(input_data) if isinstance(input_data, str) else input_data
            output_name = kwargs.get('output_name', 'merged_metabolic_pathways')
            output_file = kwargs.get('output_file', input_dir / f"{output_name}.fasta")
            
            self.logger.info(f"Starting FASTA file merging", 
                           input_dir=str(input_dir),
                           output_file=str(output_file))
            
            # Step 1: Find FASTA files
            fasta_files = self._find_fasta_files(input_dir)
            self.merge_stats["files_found"] = len(fasta_files)
            
            if not fasta_files:
                result.add_error("No FASTA files found in directory")
                return result
            
            # Step 2: Merge files
            merge_result = self._merge_fasta_files(fasta_files, Path(output_file))
            result = result.merge(merge_result)
            
            # Update metadata
            result.metadata.update({
                "merge_stats": self.merge_stats,
                "input_directory": str(input_dir),
                "output_file": str(output_file),
                "deduplication_stats": {
                    "duplicates_removed": self.deduplicator.get_duplicate_count(),
                    "deduplication_enabled": self.merge_config.deduplicate_sequences,
                    "deduplication_method": self.merge_config.deduplicate_by
                }
            })
            
            self.logger.info(f"FASTA merging completed",
                           files_processed=self.merge_stats["files_processed"],
                           sequences_merged=self.merge_stats["total_sequences_merged"])
            
        except Exception as e:
            result.add_error(f"FASTA merging failed: {str(e)}")
        
        return result
    
    def _find_fasta_files(self, input_dir: Path) -> List[FASTAFileInfo]:
        """Find and analyze FASTA files in directory."""
        fasta_files = []
        
        try:
            fasta_extensions = ['.fasta', '.fa', '.fas', '.fna', '.ffn', '.faa']
            
            for root, dirs, files in os.walk(input_dir):
                for file_name in files:
                    file_path = Path(root) / file_name
                    
                    if any(file_name.lower().endswith(ext) for ext in fasta_extensions):
                        if not file_path.is_file():
                            continue
                        
                        # Create file info
                        source_name = self.merge_config.get_source_suffix(file_path)
                        file_info = FASTAFileInfo(
                            file_path=file_path,
                            source_name=source_name
                        )
                        
                        # Skip very large files if not configured for chunk processing
                        if (not self.merge_config.chunk_processing and 
                            file_info.file_size_mb > self.merge_config.max_file_size_mb):
                            self.logger.warning(f"Skipping large file", 
                                              file=str(file_path),
                                              size_mb=file_info.file_size_mb,
                                              max_size=self.merge_config.max_file_size_mb)
                            continue
                        
                        fasta_files.append(file_info)
                        self.merge_stats["total_size_mb"] += file_info.file_size_mb
                        
                        self.logger.debug(f"Found FASTA file", 
                                        file=str(file_path),
                                        source=source_name,
                                        size_mb=file_info.file_size_mb)
            
            self.logger.info(f"Found FASTA files", 
                           count=len(fasta_files),
                           total_size_mb=self.merge_stats["total_size_mb"])
            
            return fasta_files
            
        except Exception as e:
            self.logger.error(f"Failed to find FASTA files", error=str(e))
            raise
    
    def _merge_fasta_files(self, fasta_files: List[FASTAFileInfo], output_file: Path) -> ProcessingResult:
        """Merge FASTA files into single output file."""
        result = ProcessingResult(success=True)
        
        try:
            output_file.parent.mkdir(parents=True, exist_ok=True)
            
            with open(output_file, 'w') as outfile:
                for file_info in fasta_files:
                    file_result = self._process_single_fasta(file_info, outfile)
                    
                    if file_result.success:
                        self.merge_stats["files_processed"] += 1
                        
                        # Update file statistics from processing result
                        file_stats = file_result.metadata
                        file_info.sequence_count = file_stats.get("sequences_processed", 0)
                        file_info.valid_sequences = file_stats.get("valid_sequences", 0)
                        file_info.invalid_sequences = file_stats.get("invalid_sequences", 0)
                        
                        self.merge_stats["total_sequences_processed"] += file_info.sequence_count
                        self.merge_stats["total_sequences_merged"] += file_info.valid_sequences
                        self.merge_stats["sequences_failed_validation"] += file_info.invalid_sequences
                    else:
                        self.merge_stats["files_failed"] += 1
                        result.errors.extend(file_result.errors)
                        result.warnings.extend(file_result.warnings)
            
            # Add deduplication statistics
            self.merge_stats["sequences_deduplicated"] = self.deduplicator.get_duplicate_count()
            
            result.metadata.update({
                "output_file": str(output_file),
                "files_merged": self.merge_stats["files_processed"],
                "total_sequences": self.merge_stats["total_sequences_merged"],
                "duplicates_removed": self.merge_stats["sequences_deduplicated"]
            })
            
            self.logger.info(f"FASTA files merged successfully", 
                           output=str(output_file),
                           files=self.merge_stats["files_processed"],
                           sequences=self.merge_stats["total_sequences_merged"])
            
        except Exception as e:
            result.add_error(f"Failed to merge FASTA files: {str(e)}")
        
        return result
    
    def _process_single_fasta(self, file_info: FASTAFileInfo, output_handle) -> ProcessingResult:
        """Process single FASTA file and write to output."""
        result = ProcessingResult(success=True)
        
        try:
            self.logger.debug(f"Processing FASTA file", file=str(file_info.file_path))
            
            sequences_processed = 0
            valid_sequences = 0
            invalid_sequences = 0
            
            for sequence_record in SeqIO.parse(file_info.file_path, "fasta"):
                sequences_processed += 1
                
                # Validate sequence
                is_valid, validation_errors = self.sequence_validator.validate_sequence(sequence_record)
                
                if not is_valid:
                    invalid_sequences += 1
                    for error in validation_errors:
                        self.logger.warning(error, file=str(file_info.file_path))
                    continue
                
                # Check for duplicates
                if self.deduplicator.is_duplicate(sequence_record):
                    continue
                
                # Modify header if configured
                if self.merge_config.add_source_suffix:
                    new_header = f"{sequence_record.description}{self.merge_config.source_separator}{file_info.source_name}"
                    sequence_record.description = new_header
                    sequence_record.id = sequence_record.id  # Keep original ID
                
                # Write to output
                SeqIO.write(sequence_record, output_handle, "fasta")
                valid_sequences += 1
            
            result.metadata.update({
                "sequences_processed": sequences_processed,
                "valid_sequences": valid_sequences,
                "invalid_sequences": invalid_sequences,
                "file_path": str(file_info.file_path)
            })
            
            self.logger.debug(f"Processed FASTA file successfully", 
                            file=str(file_info.file_path),
                            sequences=valid_sequences,
                            invalid=invalid_sequences)
            
        except Exception as e:
            result.add_error(f"Failed to process FASTA file {file_info.file_path}: {str(e)}")
            self.logger.error(f"Failed to process FASTA file", 
                            file=str(file_info.file_path),
                            error=str(e))
        
        return result
    
    def get_merge_statistics(self) -> Dict[str, Any]:
        """Get comprehensive merge statistics."""
        base_stats = self.get_statistics()
        
        success_rate = 0.0
        if self.merge_stats["files_found"] > 0:
            success_rate = self.merge_stats["files_processed"] / self.merge_stats["files_found"]
        
        dedup_rate = 0.0
        if self.merge_stats["total_sequences_processed"] > 0:
            dedup_rate = self.merge_stats["sequences_deduplicated"] / self.merge_stats["total_sequences_processed"]
        
        return {
            **base_stats,
            **self.merge_stats,
            "success_rate": success_rate,
            "deduplication_rate": dedup_rate,
            "sequences_per_file": (
                self.merge_stats["total_sequences_merged"] / 
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
    """Main entry point for FASTA file merging."""
    # Configuration
    input_dir = Path("/groups/itay_mayrose/alongonda/datasets/MIBIG/plant_mgcs/fasta_files")
    output_name = "merged_metabolic_pathways"
    
    # Setup configuration
    config = BioinformaticsConfig(
        input_dir=input_dir,
        output_dir=input_dir,
        max_workers=1,  # Sequential processing for file writing
        log_level="INFO"
    )
    
    merge_config = FASTAMergeConfig(
        add_source_suffix=True,
        source_separator="$",
        deduplicate_sequences=True,
        deduplicate_by="sequence",
        validate_sequences=True
    )
    
    # Initialize merger
    try:
        merger = FASTAFileMerger(
            config=config,
            merge_config=merge_config
        )
        
        # Run merging
        result = merger.run(input_dir, output_name=output_name)
        
        # Display results
        if result.success:
            print("✅ FASTA file merging completed successfully!")
            stats = merger.get_merge_statistics()
            print(f"   Files processed: {stats['files_processed']}")
            print(f"   Sequences merged: {stats['total_sequences_merged']}")
            print(f"   Duplicates removed: {stats['sequences_deduplicated']}")
            print(f"   Success rate: {stats['success_rate']:.2%}")
            print(f"   Processing time: {result.processing_time:.2f} seconds")
            print(f"   Output file: {result.metadata.get('output_file')}")
        else:
            print("❌ FASTA file merging failed!")
            for error in result.errors:
                print(f"   Error: {error}")
        
        # Display warnings
        for warning in result.warnings:
            print(f"   Warning: {warning}")
        
    except Exception as e:
        print(f"❌ Critical error: {e}")
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())