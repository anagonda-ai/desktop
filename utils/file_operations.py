"""
File operation utilities for bioinformatics data handling.

This module provides standardized file I/O operations for common bioinformatics
file formats including FASTA, CSV, GFF, and more.
"""

import csv
import json
import gzip
import pickle
from typing import Any, Dict, List, Optional, Union, Iterator, Tuple
from pathlib import Path
from dataclasses import dataclass
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging

from ..core.types import SequenceData, SequenceType, FileFormat, FilePath
from ..core.exceptions import FileOperationError, ValidationError, ErrorCode
from ..core.base import FileProcessor


@dataclass
class FileMetadata:
    """Metadata for processed files."""
    file_path: Path
    file_format: str
    file_size: int
    record_count: int
    created_by: str
    timestamp: str
    checksum: Optional[str] = None


class FileManager:
    """
    Central file management utilities.
    
    Provides consistent file operations with proper error handling,
    validation, and metadata tracking.
    """
    
    def __init__(self, base_dir: Optional[FilePath] = None):
        self.base_dir = Path(base_dir) if base_dir else Path.cwd()
        self.logger = logging.getLogger(__name__)
    
    def ensure_directory(self, path: FilePath) -> Path:
        """Ensure directory exists, create if necessary."""
        dir_path = Path(path)
        dir_path.mkdir(parents=True, exist_ok=True)
        return dir_path
    
    def get_file_info(self, file_path: FilePath) -> FileMetadata:
        """Get comprehensive file information."""
        path = Path(file_path)
        
        if not path.exists():
            raise FileOperationError(
                f"File not found: {path}",
                file_path=str(path),
                error_code=ErrorCode.FILE_NOT_FOUND
            )
        
        # Detect format
        format_detector = FileFormatDetector()
        file_format = format_detector.detect_format(path)
        
        # Count records based on format
        record_count = self._count_records(path, file_format)
        
        return FileMetadata(
            file_path=path,
            file_format=file_format,
            file_size=path.stat().st_size,
            record_count=record_count,
            created_by=self.__class__.__name__,
            timestamp=pd.Timestamp.now().isoformat()
        )
    
    def _count_records(self, file_path: Path, file_format: str) -> int:
        """Count records in file based on format."""
        try:
            if file_format in ['fasta', 'fastq']:
                return sum(1 for _ in SeqIO.parse(file_path, file_format))
            elif file_format in ['csv', 'tsv']:
                separator = ',' if file_format == 'csv' else '\t'
                with open(file_path, 'r') as f:
                    return sum(1 for _ in csv.reader(f, delimiter=separator)) - 1  # Subtract header
            elif file_format == 'json':
                with open(file_path, 'r') as f:
                    data = json.load(f)
                    return len(data) if isinstance(data, list) else 1
            else:
                # For unknown formats, count lines
                with open(file_path, 'r') as f:
                    return sum(1 for _ in f)
        except Exception as e:
            self.logger.warning(f"Could not count records in {file_path}: {e}")
            return 0
    
    def backup_file(self, file_path: FilePath, backup_suffix: str = ".bak") -> Path:
        """Create backup of file."""
        source_path = Path(file_path)
        backup_path = source_path.with_suffix(source_path.suffix + backup_suffix)
        
        import shutil
        shutil.copy2(source_path, backup_path)
        
        self.logger.info(f"Created backup: {backup_path}")
        return backup_path
    
    def safe_write(self, file_path: FilePath, content: str, backup: bool = True) -> Path:
        """Safely write content to file with optional backup."""
        path = Path(file_path)
        
        # Create backup if file exists
        if backup and path.exists():
            self.backup_file(path)
        
        # Ensure directory exists
        self.ensure_directory(path.parent)
        
        # Write to temporary file first
        temp_path = path.with_suffix(path.suffix + '.tmp')
        try:
            with open(temp_path, 'w') as f:
                f.write(content)
            
            # Move temporary file to final location
            temp_path.replace(path)
            
            self.logger.info(f"Successfully wrote file: {path}")
            return path
            
        except Exception as e:
            # Clean up temporary file
            if temp_path.exists():
                temp_path.unlink()
            raise FileOperationError(
                f"Failed to write file: {e}",
                file_path=str(path),
                operation="write",
                original_error=e
            )


class FileFormatDetector:
    """Detect file formats based on extension and content."""
    
    def __init__(self):
        self.extension_map = {
            '.fasta': 'fasta',
            '.fa': 'fasta',
            '.fas': 'fasta',
            '.fna': 'fasta',
            '.ffn': 'fasta',
            '.faa': 'fasta',
            '.fastq': 'fastq',
            '.fq': 'fastq',
            '.csv': 'csv',
            '.tsv': 'tsv',
            '.txt': 'text',
            '.json': 'json',
            '.xml': 'xml',
            '.gff': 'gff',
            '.gff3': 'gff3',
            '.gtf': 'gtf',
            '.gb': 'genbank',
            '.gbk': 'genbank',
            '.embl': 'embl',
            '.bed': 'bed',
            '.sam': 'sam',
            '.bam': 'bam',
            '.vcf': 'vcf'
        }
    
    def detect_format(self, file_path: FilePath) -> str:
        """Detect file format from extension and content."""
        path = Path(file_path)
        
        # Handle compressed files
        if path.suffix == '.gz':
            base_path = path.with_suffix('')
            extension = base_path.suffix
        else:
            extension = path.suffix
        
        # Try extension-based detection first
        detected_format = self.extension_map.get(extension.lower())
        
        if detected_format:
            return detected_format
        
        # Fall back to content-based detection
        return self._detect_from_content(path)
    
    def _detect_from_content(self, file_path: Path) -> str:
        """Detect format from file content."""
        try:
            # Read first few lines
            opener = gzip.open if file_path.suffix == '.gz' else open
            mode = 'rt' if file_path.suffix == '.gz' else 'r'
            
            with opener(file_path, mode) as f:
                first_lines = [f.readline().strip() for _ in range(5)]
            
            # FASTA format detection
            if any(line.startswith('>') for line in first_lines):
                return 'fasta'
            
            # FASTQ format detection
            if first_lines[0].startswith('@') and first_lines[2].startswith('+'):
                return 'fastq'
            
            # GFF/GTF detection
            if any('\t' in line and len(line.split('\t')) >= 8 for line in first_lines):
                return 'gff'
            
            # JSON detection
            if first_lines[0].startswith('{') or first_lines[0].startswith('['):
                return 'json'
            
            # CSV detection (comma-separated)
            if any(',' in line for line in first_lines):
                return 'csv'
            
            # TSV detection (tab-separated)
            if any('\t' in line for line in first_lines):
                return 'tsv'
            
            return 'text'  # Default fallback
            
        except Exception:
            return 'unknown'


class CSVHandler:
    """Standardized CSV file operations."""
    
    def __init__(self, encoding: str = 'utf-8'):
        self.encoding = encoding
        self.logger = logging.getLogger(__name__)
    
    def read_csv(
        self, 
        file_path: FilePath, 
        **kwargs
    ) -> pd.DataFrame:
        """Read CSV file with error handling."""
        try:
            df = pd.read_csv(file_path, encoding=self.encoding, **kwargs)
            self.logger.debug(f"Read CSV with {len(df)} rows from {file_path}")
            return df
        except Exception as e:
            raise FileOperationError(
                f"Failed to read CSV file: {e}",
                file_path=str(file_path),
                operation="read_csv",
                original_error=e
            )
    
    def write_csv(
        self, 
        df: pd.DataFrame, 
        file_path: FilePath, 
        index: bool = False,
        **kwargs
    ) -> Path:
        """Write DataFrame to CSV with error handling."""
        path = Path(file_path)
        try:
            # Ensure directory exists
            path.parent.mkdir(parents=True, exist_ok=True)
            
            df.to_csv(path, index=index, encoding=self.encoding, **kwargs)
            self.logger.debug(f"Wrote CSV with {len(df)} rows to {path}")
            return path
        except Exception as e:
            raise FileOperationError(
                f"Failed to write CSV file: {e}",
                file_path=str(path),
                operation="write_csv",
                original_error=e
            )
    
    def append_csv(
        self, 
        df: pd.DataFrame, 
        file_path: FilePath, 
        **kwargs
    ) -> Path:
        """Append DataFrame to existing CSV file."""
        path = Path(file_path)
        
        # Check if file exists and has headers
        write_header = not path.exists()
        
        try:
            df.to_csv(
                path, 
                mode='a', 
                header=write_header, 
                index=False,
                encoding=self.encoding,
                **kwargs
            )
            self.logger.debug(f"Appended {len(df)} rows to {path}")
            return path
        except Exception as e:
            raise FileOperationError(
                f"Failed to append to CSV file: {e}",
                file_path=str(path),
                operation="append_csv",
                original_error=e
            )


class FASTAHandler:
    """Standardized FASTA file operations."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def read_fasta(
        self, 
        file_path: FilePath, 
        sequence_type: Optional[SequenceType] = None
    ) -> List[SequenceData]:
        """Read FASTA file and return sequence data objects."""
        path = Path(file_path)
        sequences = []
        
        try:
            for record in SeqIO.parse(path, "fasta"):
                seq_type = sequence_type or self._detect_sequence_type(str(record.seq))
                
                seq_data = SequenceData(
                    sequence_id=record.id,
                    sequence=str(record.seq),
                    sequence_type=seq_type,
                    description=record.description,
                    organism=self._extract_organism(record.description)
                )
                sequences.append(seq_data)
            
            self.logger.debug(f"Read {len(sequences)} sequences from {path}")
            return sequences
            
        except Exception as e:
            raise FileOperationError(
                f"Failed to read FASTA file: {e}",
                file_path=str(path),
                operation="read_fasta",
                original_error=e
            )
    
    def write_fasta(
        self, 
        sequences: List[SequenceData], 
        file_path: FilePath,
        line_length: int = 80
    ) -> Path:
        """Write sequence data to FASTA file."""
        path = Path(file_path)
        
        try:
            # Ensure directory exists
            path.parent.mkdir(parents=True, exist_ok=True)
            
            records = []
            for seq_data in sequences:
                record = SeqRecord(
                    Seq(seq_data.sequence),
                    id=seq_data.sequence_id,
                    description=seq_data.description or ""
                )
                records.append(record)
            
            SeqIO.write(records, path, "fasta")
            self.logger.debug(f"Wrote {len(sequences)} sequences to {path}")
            return path
            
        except Exception as e:
            raise FileOperationError(
                f"Failed to write FASTA file: {e}",
                file_path=str(path),
                operation="write_fasta",
                original_error=e
            )
    
    def _detect_sequence_type(self, sequence: str) -> SequenceType:
        """Detect sequence type from sequence content."""
        sequence = sequence.upper().replace('-', '').replace('N', '')
        
        if not sequence:
            return SequenceType.DNA
        
        # Count nucleotides
        nucleotides = set('ATCG')
        protein_chars = set('ARNDCEQGHILKMFPSTWYV')
        
        nucleotide_count = sum(1 for char in sequence if char in nucleotides)
        nucleotide_ratio = nucleotide_count / len(sequence)
        
        # If >90% are standard nucleotides, it's likely DNA/RNA
        if nucleotide_ratio > 0.9:
            if 'U' in sequence:
                return SequenceType.RNA
            else:
                return SequenceType.DNA
        else:
            return SequenceType.PROTEIN
    
    def _extract_organism(self, description: str) -> Optional[str]:
        """Extract organism name from FASTA description."""
        if not description:
            return None
        
        # Common patterns for organism names in FASTA headers
        patterns = [
            r'\[([^\]]+)\]',  # [Organism name]
            r'OS=([^=]+?)(?:\s+[A-Z]{2}=|$)',  # OS=Organism name
        ]
        
        import re
        for pattern in patterns:
            match = re.search(pattern, description)
            if match:
                return match.group(1).strip()
        
        return None
    
    def merge_fasta_files(
        self, 
        input_files: List[FilePath], 
        output_file: FilePath,
        remove_duplicates: bool = True
    ) -> Path:
        """Merge multiple FASTA files into one."""
        all_sequences = []
        seen_ids = set()
        
        for file_path in input_files:
            sequences = self.read_fasta(file_path)
            
            for seq in sequences:
                if remove_duplicates:
                    if seq.sequence_id not in seen_ids:
                        all_sequences.append(seq)
                        seen_ids.add(seq.sequence_id)
                else:
                    all_sequences.append(seq)
        
        return self.write_fasta(all_sequences, output_file)


class SequenceFileReader(FileProcessor):
    """Unified reader for various sequence file formats."""
    
    def __init__(self):
        super().__init__()
        self.supported_formats = ['fasta', 'fastq', 'genbank', 'embl']
        self.fasta_handler = FASTAHandler()
    
    def validate_input(self, file_path: FilePath) -> None:
        """Validate input file."""
        path = self.validate_file_exists(file_path, "sequence file")
        file_format = self.detect_file_format(path)
        self.validate_file_format(path, file_format)
    
    def process(self, file_path: FilePath, **kwargs) -> List[SequenceData]:
        """Read sequences from file."""
        path = Path(file_path)
        file_format = self.detect_file_format(path)
        
        if file_format == 'fasta':
            return self.fasta_handler.read_fasta(path, **kwargs)
        elif file_format == 'fastq':
            return self._read_fastq(path, **kwargs)
        elif file_format in ['genbank', 'gb']:
            return self._read_genbank(path, **kwargs)
        else:
            raise ValidationError(
                f"Unsupported sequence format: {file_format}",
                error_code=ErrorCode.FILE_FORMAT_INVALID
            )
    
    def _read_fastq(self, file_path: Path, **kwargs) -> List[SequenceData]:
        """Read FASTQ file."""
        sequences = []
        
        for record in SeqIO.parse(file_path, "fastq"):
            seq_data = SequenceData(
                sequence_id=record.id,
                sequence=str(record.seq),
                sequence_type=SequenceType.DNA,  # FASTQ is typically DNA
                description=record.description
            )
            sequences.append(seq_data)
        
        return sequences
    
    def _read_genbank(self, file_path: Path, **kwargs) -> List[SequenceData]:
        """Read GenBank file."""
        sequences = []
        
        for record in SeqIO.parse(file_path, "genbank"):
            seq_data = SequenceData(
                sequence_id=record.id,
                sequence=str(record.seq),
                sequence_type=SequenceType.DNA,  # GenBank is typically DNA
                description=record.description,
                organism=record.annotations.get('organism')
            )
            sequences.append(seq_data)
        
        return sequences


class SequenceFileWriter(FileProcessor):
    """Unified writer for various sequence file formats."""
    
    def __init__(self):
        super().__init__()
        self.supported_formats = ['fasta', 'fastq', 'genbank']
        self.fasta_handler = FASTAHandler()
    
    def validate_input(self, data: List[SequenceData]) -> None:
        """Validate sequence data."""
        if not isinstance(data, list):
            raise ValidationError("Input must be a list of SequenceData objects")
        
        if not data:
            raise ValidationError("Input list cannot be empty")
        
        for i, seq in enumerate(data):
            if not isinstance(seq, SequenceData):
                raise ValidationError(f"Item {i} is not a SequenceData object")
    
    def process(
        self, 
        sequences: List[SequenceData], 
        output_file: FilePath,
        output_format: str = 'fasta',
        **kwargs
    ) -> Path:
        """Write sequences to file."""
        if output_format == 'fasta':
            return self.fasta_handler.write_fasta(sequences, output_file, **kwargs)
        else:
            raise ValidationError(
                f"Unsupported output format: {output_format}",
                error_code=ErrorCode.FILE_FORMAT_INVALID
            )