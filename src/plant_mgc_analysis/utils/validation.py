"""
Comprehensive Data Validation Layer for Plant MGC Analysis Pipeline.

This module provides robust validation for all data types, analysis parameters,
and input files used throughout the bioinformatics pipeline with enhanced
object-oriented design and comprehensive error handling.
"""

import re
import os
from typing import Dict, List, Optional, Any, Union, Set, Tuple, Type
from dataclasses import dataclass, field
from pathlib import Path
from abc import ABC, abstractmethod
import json

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from loguru import logger

from ..core.base import LoggerMixin
from ..core.types import AnalysisType, DatabaseType, GeneInfo, GenomeData, MGCCandidate, PathLike
from ..core.exceptions import ValidationError, ConfigurationError, FileSystemError


@dataclass
class ValidationRule:
    """Represents a validation rule with metadata."""
    
    name: str
    description: str
    error_message: str
    is_critical: bool = True
    
    def __post_init__(self):
        """Validate rule definition."""
        if not self.name or not self.description:
            raise ConfigurationError("Validation rule must have name and description")


@dataclass
class ValidationResult:
    """Result of validation operation."""
    
    is_valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    field_name: Optional[str] = None
    field_value: Any = None
    
    def add_error(self, message: str, is_critical: bool = True) -> None:
        """Add validation error."""
        if is_critical:
            self.errors.append(message)
            self.is_valid = False
        else:
            self.warnings.append(message)
    
    def merge(self, other: 'ValidationResult') -> 'ValidationResult':
        """Merge two validation results."""
        return ValidationResult(
            is_valid=self.is_valid and other.is_valid,
            errors=self.errors + other.errors,
            warnings=self.warnings + other.warnings,
        )


class BaseValidator(ABC, LoggerMixin):
    """Abstract base class for all validators."""
    
    def __init__(self, rules: Optional[List[ValidationRule]] = None):
        """Initialize validator with optional custom rules."""
        self.rules = rules or []
        self._setup_default_rules()
    
    @abstractmethod
    def _setup_default_rules(self) -> None:
        """Setup default validation rules."""
        pass
    
    @abstractmethod
    def validate(self, data: Any, **kwargs) -> ValidationResult:
        """Validate data according to rules."""
        pass
    
    def add_rule(self, rule: ValidationRule) -> None:
        """Add custom validation rule."""
        self.rules.append(rule)
    
    def _create_error(self, message: str, field_name: str, field_value: Any) -> ValidationError:
        """Create validation error with context."""
        return ValidationError(
            message,
            field_name=field_name,
            field_value=field_value
        )


class EnhancedSequenceValidator(BaseValidator):
    """Enhanced validator for biological sequences."""
    
    def _setup_default_rules(self) -> None:
        """Setup default sequence validation rules."""
        self.rules = [
            ValidationRule(
                name="non_empty_sequence",
                description="Sequence must not be empty",
                error_message="Sequence cannot be empty"
            ),
            ValidationRule(
                name="valid_nucleotide_sequence",
                description="Nucleotide sequence must contain only valid bases",
                error_message="Invalid nucleotide characters found"
            ),
            ValidationRule(
                name="valid_protein_sequence",
                description="Protein sequence must contain only valid amino acids",
                error_message="Invalid amino acid characters found"
            ),
            ValidationRule(
                name="minimum_length",
                description="Sequence must meet minimum length requirement",
                error_message="Sequence too short"
            ),
        ]
    
    def validate(self, data: Union[str, Seq], sequence_type: str = "auto", **kwargs) -> ValidationResult:
        """Enhanced sequence validation with auto-detection."""
        result = ValidationResult(is_valid=True)
        
        # Convert to string if needed
        if hasattr(data, 'data'):
            sequence = str(data.data)
        elif hasattr(data, '__str__'):
            sequence = str(data)
        else:
            result.add_error("Invalid sequence type", True)
            return result
        
        sequence = sequence.upper().strip()
        
        # Check if empty
        if not sequence:
            result.add_error("Sequence cannot be empty", True)
            return result
        
        # Auto-detect sequence type if needed
        if sequence_type == "auto":
            sequence_type = self._detect_sequence_type(sequence)
        
        # Validate based on sequence type
        if sequence_type == "nucleotide":
            self._validate_nucleotide_sequence(sequence, result, **kwargs)
        elif sequence_type == "protein":
            self._validate_protein_sequence(sequence, result, **kwargs)
        else:
            result.add_error(f"Unknown sequence type: {sequence_type}", True)
        
        # Check minimum length
        min_length = kwargs.get("min_length", 1)
        if len(sequence) < min_length:
            result.add_error(f"Sequence too short (minimum {min_length} characters)", True)
        
        # Check maximum length
        max_length = kwargs.get("max_length")
        if max_length and len(sequence) > max_length:
            result.add_error(f"Sequence too long (maximum {max_length} characters)", False)
        
        return result
    
    def _detect_sequence_type(self, sequence: str) -> str:
        """Auto-detect sequence type based on content."""
        nucleotide_chars = set("ATCGNU")
        protein_chars = set("ACDEFGHIKLMNPQRSTVWY")
        
        sequence_chars = set(sequence.upper())
        
        if sequence_chars.issubset(nucleotide_chars):
            return "nucleotide"
        elif sequence_chars.issubset(protein_chars):
            return "protein"
        else:
            return "unknown"
    
    def _validate_nucleotide_sequence(self, sequence: str, result: ValidationResult, **kwargs) -> None:
        """Validate nucleotide sequence."""
        valid_chars = set("ATCGNU-")
        invalid_chars = set(sequence) - valid_chars
        
        if invalid_chars:
            result.add_error(f"Invalid nucleotide characters: {', '.join(invalid_chars)}", True)
        
        # Check GC content
        if "check_gc_content" in kwargs:
            gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
            min_gc = kwargs.get("min_gc_content", 0.0)
            max_gc = kwargs.get("max_gc_content", 1.0)
            
            if gc_content < min_gc or gc_content > max_gc:
                result.add_error(f"GC content {gc_content:.2%} outside expected range", False)
    
    def _validate_protein_sequence(self, sequence: str, result: ValidationResult, **kwargs) -> None:
        """Validate protein sequence."""
        valid_chars = set("ACDEFGHIKLMNPQRSTVWY*-")
        invalid_chars = set(sequence) - valid_chars
        
        if invalid_chars:
            result.add_error(f"Invalid amino acid characters: {', '.join(invalid_chars)}", True)
        
        # Check for stop codons in middle
        if "check_internal_stops" in kwargs:
            internal_stops = sequence[:-1].count('*')
            if internal_stops > 0:
                result.add_error(f"Found {internal_stops} internal stop codons", False)


class EnhancedFileValidator(BaseValidator):
    """Enhanced validator for input files with comprehensive format checking."""
    
    def _setup_default_rules(self) -> None:
        """Setup default file validation rules."""
        self.rules = [
            ValidationRule(
                name="file_exists",
                description="File must exist",
                error_message="File does not exist"
            ),
            ValidationRule(
                name="file_readable",
                description="File must be readable",
                error_message="File is not readable"
            ),
            ValidationRule(
                name="file_not_empty",
                description="File must not be empty",
                error_message="File is empty"
            ),
            ValidationRule(
                name="valid_format",
                description="File must be in expected format",
                error_message="Invalid file format"
            ),
        ]
    
    def validate(self, data: Union[str, Path], file_type: str = "auto", **kwargs) -> ValidationResult:
        """Enhanced file validation with auto-detection."""
        result = ValidationResult(is_valid=True)
        
        # Convert to Path object
        file_path = Path(data)
        
        # Check if file exists
        if not file_path.exists():
            result.add_error(f"File does not exist: {file_path}", True)
            return result
        
        # Check if file is readable
        if not file_path.is_file():
            result.add_error(f"Path is not a file: {file_path}", True)
            return result
        
        # Check file permissions
        if not os.access(file_path, os.R_OK):
            result.add_error(f"File is not readable: {file_path}", True)
            return result
        
        # Check if file is empty
        if file_path.stat().st_size == 0:
            result.add_error(f"File is empty: {file_path}", True)
            return result
        
        # Auto-detect file type if needed
        if file_type == "auto":
            file_type = self._detect_file_type(file_path)
        
        # Validate based on file type
        if file_type == "fasta":
            self._validate_fasta_file(file_path, result, **kwargs)
        elif file_type == "gff":
            self._validate_gff_file(file_path, result, **kwargs)
        elif file_type == "csv":
            self._validate_csv_file(file_path, result, **kwargs)
        elif file_type == "json":
            self._validate_json_file(file_path, result, **kwargs)
        
        return result
    
    def _detect_file_type(self, file_path: Path) -> str:
        """Auto-detect file type based on extension and content."""
        extension = file_path.suffix.lower()
        
        extension_map = {
            '.fa': 'fasta', '.fas': 'fasta', '.fasta': 'fasta',
            '.gff': 'gff', '.gff3': 'gff',
            '.csv': 'csv', '.tsv': 'csv',
            '.json': 'json',
        }
        
        return extension_map.get(extension, "unknown")
    
    def _validate_fasta_file(self, file_path: Path, result: ValidationResult, **kwargs) -> None:
        """Enhanced FASTA validation."""
        try:
            sequences = list(SeqIO.parse(file_path, "fasta"))
            
            if not sequences:
                result.add_error("No sequences found in FASTA file", True)
                return
            
            # Check sequence count
            min_sequences = kwargs.get("min_sequences", 1)
            if len(sequences) < min_sequences:
                result.add_error(f"Too few sequences (minimum {min_sequences})", True)
            
            # Validate individual sequences
            sequence_validator = EnhancedSequenceValidator()
            
            for i, seq_record in enumerate(sequences):
                if not seq_record.id:
                    result.add_error(f"Sequence {i+1} missing identifier", True)
                
                if not str(seq_record.seq):
                    result.add_error(f"Sequence {i+1} is empty", True)
                
                # Validate sequence content
                seq_result = sequence_validator.validate(seq_record.seq, **kwargs)
                if not seq_result.is_valid:
                    for error in seq_result.errors:
                        result.add_error(f"Sequence {i+1}: {error}", True)
        
        except Exception as e:
            result.add_error(f"Error parsing FASTA file: {e}", True)
    
    def _validate_gff_file(self, file_path: Path, result: ValidationResult, **kwargs) -> None:
        """Enhanced GFF validation."""
        try:
            with open(file_path, 'r') as f:
                line_count = 0
                feature_count = 0
                
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    line_count += 1
                    
                    # Skip comments and empty lines
                    if line.startswith('#') or not line:
                        continue
                    
                    # Parse GFF line
                    fields = line.split('\t')
                    if len(fields) != 9:
                        result.add_error(f"Line {line_num}: Invalid number of fields (expected 9)", True)
                        continue
                    
                    feature_count += 1
                    
                    # Validate coordinates
                    try:
                        start = int(fields[3])
                        end = int(fields[4])
                        
                        if start > end:
                            result.add_error(f"Line {line_num}: Start > End", True)
                        if start < 1:
                            result.add_error(f"Line {line_num}: Invalid start coordinate", True)
                    
                    except ValueError:
                        result.add_error(f"Line {line_num}: Invalid coordinates", True)
                
                if feature_count == 0:
                    result.add_error("No features found in GFF file", True)
        
        except Exception as e:
            result.add_error(f"Error parsing GFF file: {e}", True)
    
    def _validate_csv_file(self, file_path: Path, result: ValidationResult, **kwargs) -> None:
        """Enhanced CSV validation."""
        try:
            df = pd.read_csv(file_path)
            
            if df.empty:
                result.add_error("CSV file is empty", True)
                return
            
            # Check required columns
            required_columns = kwargs.get("required_columns", [])
            missing_columns = set(required_columns) - set(df.columns)
            
            if missing_columns:
                result.add_error(f"Missing required columns: {', '.join(missing_columns)}", True)
            
            # Check for duplicate columns
            if len(df.columns) != len(set(df.columns)):
                result.add_error("Duplicate column names found", True)
        
        except Exception as e:
            result.add_error(f"Error parsing CSV file: {e}", True)
    
    def _validate_json_file(self, file_path: Path, result: ValidationResult, **kwargs) -> None:
        """Enhanced JSON validation."""
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            # Check required keys
            required_keys = kwargs.get("required_keys", [])
            if isinstance(data, dict) and required_keys:
                missing_keys = set(required_keys) - set(data.keys())
                
                if missing_keys:
                    result.add_error(f"Missing required keys: {', '.join(missing_keys)}", True)
        
        except json.JSONDecodeError as e:
            result.add_error(f"Invalid JSON format: {e}", True)
        except Exception as e:
            result.add_error(f"Error parsing JSON file: {e}", True)


class CompositeValidator:
    """Composite validator that combines multiple validators for complete validation."""
    
    def __init__(self):
        """Initialize composite validator."""
        self.sequence_validator = EnhancedSequenceValidator()
        self.file_validator = EnhancedFileValidator()
    
    def validate_genome_data(self, genome_data: GenomeData) -> ValidationResult:
        """Validate complete genome data."""
        result = ValidationResult(is_valid=True)
        
        # Validate genome file
        if genome_data.genome_file:
            file_result = self.file_validator.validate(
                genome_data.genome_file,
                file_type="fasta",
                min_sequences=1
            )
            result = result.merge(file_result)
        
        # Validate protein file
        if genome_data.protein_file:
            file_result = self.file_validator.validate(
                genome_data.protein_file,
                file_type="fasta",
                min_sequences=1,
                sequence_type="protein"
            )
            result = result.merge(file_result)
        
        # Validate annotation file
        if genome_data.annotation_file:
            file_result = self.file_validator.validate(
                genome_data.annotation_file,
                file_type="gff"
            )
            result = result.merge(file_result)
        
        return result


def validate_input_files(file_paths: List[PathLike]) -> None:
    """
    Validate that input files exist and are readable.
    
    Args:
        file_paths: List of file paths to validate
        
    Raises:
        ValidationError: If any file is invalid
        FileSystemError: If file system operations fail
    """
    for file_path in file_paths:
        path = Path(file_path)
        
        if not path.exists():
            raise FileSystemError(
                f"File does not exist: {path}",
                file_path=str(path),
                operation="validate_existence"
            )
        
        if not path.is_file():
            raise ValidationError(
                f"Path is not a file: {path}",
                field_name="file_path",
                field_value=str(path)
            )
        
        if not path.stat().st_size > 0:
            raise ValidationError(
                f"File is empty: {path}",
                field_name="file_path",
                field_value=str(path)
            )
        
        try:
            with open(path, 'r') as f:
                f.read(1)  # Try to read first character
        except PermissionError:
            raise FileSystemError(
                f"Permission denied reading file: {path}",
                file_path=str(path),
                operation="validate_readable"
            )
        except Exception as e:
            raise FileSystemError(
                f"Error reading file {path}: {e}",
                file_path=str(path),
                operation="validate_readable"
            )


def validate_fasta_file(file_path: PathLike) -> bool:
    """
    Validate FASTA file format.
    
    Args:
        file_path: Path to FASTA file
        
    Returns:
        True if valid FASTA file
        
    Raises:
        ValidationError: If file is not valid FASTA
    """
    path = Path(file_path)
    
    try:
        sequences = list(SeqIO.parse(path, "fasta"))
        
        if not sequences:
            raise ValidationError(
                f"No sequences found in FASTA file: {path}",
                field_name="fasta_file",
                field_value=str(path)
            )
        
        # Validate sequence IDs
        for i, seq in enumerate(sequences):
            if not seq.id:
                raise ValidationError(
                    f"Sequence {i+1} has no ID in FASTA file: {path}",
                    field_name="sequence_id",
                    field_value=f"sequence_{i+1}"
                )
            
            if len(seq.seq) == 0:
                raise ValidationError(
                    f"Sequence {seq.id} is empty in FASTA file: {path}",
                    field_name="sequence_length",
                    field_value=seq.id
                )
        
        logger.info(f"FASTA file validated: {path} ({len(sequences)} sequences)")
        return True
        
    except Exception as e:
        if isinstance(e, ValidationError):
            raise
        raise ValidationError(
            f"Error parsing FASTA file {path}: {e}",
            field_name="fasta_file",
            field_value=str(path)
        )


def validate_gff_file(file_path: PathLike) -> bool:
    """
    Validate GFF/GTF file format.
    
    Args:
        file_path: Path to GFF/GTF file
        
    Returns:
        True if valid GFF file
        
    Raises:
        ValidationError: If file is not valid GFF
    """
    path = Path(file_path)
    
    try:
        with open(path, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            raise ValidationError(
                f"GFF file is empty: {path}",
                field_name="gff_file",
                field_value=str(path)
            )
        
        # Check for GFF header
        if not any(line.startswith("##gff-version") for line in lines[:10]):
            logger.warning(f"GFF file missing version header: {path}")
        
        # Validate GFF format
        feature_count = 0
        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            
            # Skip comments and empty lines
            if not line or line.startswith("#"):
                continue
            
            parts = line.split('\t')
            if len(parts) != 9:
                raise ValidationError(
                    f"Invalid GFF format at line {line_num}: expected 9 columns, got {len(parts)}",
                    field_name="gff_format",
                    field_value=f"line_{line_num}"
                )
            
            # Validate coordinates
            try:
                start = int(parts[3])
                end = int(parts[4])
                
                if start < 1 or end < 1:
                    raise ValidationError(
                        f"Invalid coordinates at line {line_num}: start and end must be >= 1",
                        field_name="coordinates",
                        field_value=f"start={start}, end={end}"
                    )
                
                if start > end:
                    raise ValidationError(
                        f"Invalid coordinates at line {line_num}: start > end",
                        field_name="coordinates",
                        field_value=f"start={start}, end={end}"
                    )
                
            except ValueError:
                raise ValidationError(
                    f"Invalid coordinate format at line {line_num}",
                    field_name="coordinates",
                    field_value=f"start={parts[3]}, end={parts[4]}"
                )
            
            feature_count += 1
        
        if feature_count == 0:
            raise ValidationError(
                f"No features found in GFF file: {path}",
                field_name="gff_features",
                field_value=str(path)
            )
        
        logger.info(f"GFF file validated: {path} ({feature_count} features)")
        return True
        
    except Exception as e:
        if isinstance(e, ValidationError):
            raise
        raise ValidationError(
            f"Error parsing GFF file {path}: {e}",
            field_name="gff_file",
            field_value=str(path)
        )


def validate_organism_name(name: str) -> bool:
    """
    Validate organism name format.
    
    Args:
        name: Organism name to validate
        
    Returns:
        True if valid organism name
        
    Raises:
        ValidationError: If name is invalid
    """
    if not name or not name.strip():
        raise ValidationError(
            "Organism name cannot be empty",
            field_name="organism_name",
            field_value=name
        )
    
    # Check for reasonable length
    if len(name) > 100:
        raise ValidationError(
            "Organism name too long (max 100 characters)",
            field_name="organism_name",
            field_value=name
        )
    
    # Check for invalid characters
    if re.search(r'[<>:"/\\|?*]', name):
        raise ValidationError(
            "Organism name contains invalid characters",
            field_name="organism_name",
            field_value=name
        )
    
    return True


def validate_analysis_parameters(
    analysis_type: str, 
    parameters: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Validate analysis parameters for specific analysis type.
    
    Args:
        analysis_type: Type of analysis
        parameters: Parameters to validate
        
    Returns:
        Validated parameters dictionary
        
    Raises:
        ValidationError: If parameters are invalid
    """
    validated = {}
    
    if analysis_type == "sliding_window":
        # Validate window size
        window_size = parameters.get("window_size", 50000)
        if not isinstance(window_size, int) or window_size <= 0:
            raise ValidationError(
                "Window size must be a positive integer",
                field_name="window_size",
                field_value=window_size,
                expected_type=int
            )
        validated["window_size"] = window_size
        
        # Validate step size
        step_size = parameters.get("step_size", 10000)
        if not isinstance(step_size, int) or step_size <= 0:
            raise ValidationError(
                "Step size must be a positive integer",
                field_name="step_size",
                field_value=step_size,
                expected_type=int
            )
        validated["step_size"] = step_size
        
        # Validate min genes
        min_genes = parameters.get("min_genes", 3)
        if not isinstance(min_genes, int) or min_genes < 1:
            raise ValidationError(
                "Minimum genes must be a positive integer",
                field_name="min_genes",
                field_value=min_genes,
                expected_type=int
            )
        validated["min_genes"] = min_genes
    
    elif analysis_type == "blast_search":
        # Validate e-value
        evalue = parameters.get("evalue", 1e-5)
        if not isinstance(evalue, (int, float)) or evalue <= 0:
            raise ValidationError(
                "E-value must be a positive number",
                field_name="evalue",
                field_value=evalue,
                expected_type=float
            )
        validated["evalue"] = float(evalue)
        
        # Validate max targets
        max_targets = parameters.get("max_targets", 100)
        if not isinstance(max_targets, int) or max_targets < 1:
            raise ValidationError(
                "Max targets must be a positive integer",
                field_name="max_targets",
                field_value=max_targets,
                expected_type=int
            )
        validated["max_targets"] = max_targets
    
    # Add other analysis type validations as needed
    
    return validated


def validate_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate configuration dictionary.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        Validated configuration
        
    Raises:
        ValidationError: If configuration is invalid
    """
    validated = {}
    
    # Validate required fields
    required_fields = ["data_dir", "output_dir"]
    for field in required_fields:
        if field not in config:
            raise ValidationError(
                f"Required configuration field missing: {field}",
                field_name=field
            )
    
    # Validate directories
    for dir_field in ["data_dir", "output_dir", "cache_dir"]:
        if dir_field in config:
            dir_path = Path(config[dir_field])
            try:
                dir_path.mkdir(parents=True, exist_ok=True)
                validated[dir_field] = dir_path
            except Exception as e:
                raise ValidationError(
                    f"Cannot create directory {dir_field}: {e}",
                    field_name=dir_field,
                    field_value=str(dir_path)
                )
    
    # Validate numeric parameters
    numeric_fields = {
        "max_workers": int,
        "chunk_size": int,
        "max_memory_gb": float,
    }
    
    for field, expected_type in numeric_fields.items():
        if field in config:
            value = config[field]
            if not isinstance(value, expected_type) or value <= 0:
                raise ValidationError(
                    f"Invalid {field}: must be positive {expected_type.__name__}",
                    field_name=field,
                    field_value=value,
                    expected_type=expected_type
                )
            validated[field] = value
    
    # Copy other valid fields
    for key, value in config.items():
        if key not in validated:
            validated[key] = value
    
    return validated


def validate_sequence_ids(sequence_ids: List[str]) -> List[str]:
    """
    Validate sequence IDs.
    
    Args:
        sequence_ids: List of sequence IDs to validate
        
    Returns:
        List of validated sequence IDs
        
    Raises:
        ValidationError: If any ID is invalid
    """
    if not sequence_ids:
        raise ValidationError(
            "Sequence ID list cannot be empty",
            field_name="sequence_ids",
            field_value=sequence_ids
        )
    
    validated = []
    for i, seq_id in enumerate(sequence_ids):
        if not seq_id or not seq_id.strip():
            raise ValidationError(
                f"Sequence ID {i+1} is empty",
                field_name="sequence_id",
                field_value=seq_id
            )
        
        # Check for invalid characters
        if re.search(r'[<>:"/\\|?*\s]', seq_id):
            raise ValidationError(
                f"Sequence ID contains invalid characters: {seq_id}",
                field_name="sequence_id",
                field_value=seq_id
            )
        
        validated.append(seq_id.strip())
    
    # Check for duplicates
    if len(set(validated)) != len(validated):
        raise ValidationError(
            "Duplicate sequence IDs found",
            field_name="sequence_ids",
            field_value=sequence_ids
        )
    
    return validated