"""
File operations utilities for Plant MGC Analysis Pipeline.

This module provides high-level file operations with proper error handling,
validation, and format support for bioinformatics data.
"""

import csv
import json
import gzip
import pickle
from typing import Any, Dict, List, Optional, Union, Iterator, TextIO
from pathlib import Path
from contextlib import contextmanager
from dataclasses import asdict

import pandas as pd
from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord

from ..core.base import FileProcessor
from ..core.types import PathLike, GeneInfo, MGCCandidate
from ..core.exceptions import FileSystemError, ValidationError
from .validation import validate_fasta_file, validate_gff_file


class SafeFileOperations:
    """Safe file operations with proper error handling."""
    
    @staticmethod
    @contextmanager
    def safe_open(
        file_path: PathLike,
        mode: str = 'r',
        encoding: str = 'utf-8',
        **kwargs
    ) -> Iterator[TextIO]:
        """
        Safely open a file with proper error handling.
        
        Args:
            file_path: Path to file
            mode: File open mode
            encoding: File encoding
            **kwargs: Additional arguments for open()
            
        Yields:
            File handle
            
        Raises:
            FileSystemError: If file operations fail
        """
        path = Path(file_path)
        
        try:
            if mode.startswith('w') or mode.startswith('a'):
                # Ensure parent directory exists for write operations
                path.parent.mkdir(parents=True, exist_ok=True)
            
            # Handle compressed files
            if path.suffix == '.gz':
                file_handle = gzip.open(path, mode + 't', encoding=encoding, **kwargs)
            else:
                file_handle = open(path, mode, encoding=encoding, **kwargs)
            
            yield file_handle
            
        except PermissionError as e:
            raise FileSystemError(
                f"Permission denied accessing file: {path}",
                file_path=str(path),
                operation=f"open_{mode}"
            ) from e
        except OSError as e:
            raise FileSystemError(
                f"OS error accessing file: {path} - {e}",
                file_path=str(path),
                operation=f"open_{mode}"
            ) from e
        finally:
            if 'file_handle' in locals():
                file_handle.close()
    
    @staticmethod
    def ensure_directory(dir_path: PathLike) -> Path:
        """
        Ensure directory exists, create if necessary.
        
        Args:
            dir_path: Directory path
            
        Returns:
            Path object for directory
            
        Raises:
            FileSystemError: If directory cannot be created
        """
        path = Path(dir_path)
        
        try:
            path.mkdir(parents=True, exist_ok=True)
            return path
        except PermissionError as e:
            raise FileSystemError(
                f"Permission denied creating directory: {path}",
                file_path=str(path),
                operation="mkdir"
            ) from e
        except OSError as e:
            raise FileSystemError(
                f"OS error creating directory: {path} - {e}",
                file_path=str(path),
                operation="mkdir"
            ) from e


class FastaProcessor(FileProcessor):
    """Processor for FASTA format files."""
    
    def __init__(self, **kwargs):
        super().__init__(
            supported_formats=['.fasta', '.fa', '.fas', '.fna', '.ffn', '.faa'],
            **kwargs
        )
    
    def validate_input(self, file_path: PathLike) -> None:
        """Validate FASTA file."""
        super().validate_input(file_path)
        validate_fasta_file(file_path)
    
    def read_file(self, file_path: PathLike) -> List[SeqRecord]:
        """
        Read FASTA file and return sequences.
        
        Args:
            file_path: Path to FASTA file
            
        Returns:
            List of sequence records
        """
        self.validate_input(file_path)
        
        with SafeFileOperations.safe_open(file_path, 'r') as handle:
            sequences = list(SeqIO.parse(handle, "fasta"))
        
        self.logger.info(f"Read {len(sequences)} sequences from {file_path}")
        return sequences
    
    def write_file(self, sequences: List[SeqRecord], file_path: PathLike) -> None:
        """
        Write sequences to FASTA file.
        
        Args:
            sequences: List of sequence records
            file_path: Output file path
        """
        if not sequences:
            raise ValidationError(
                "Cannot write empty sequence list",
                field_name="sequences",
                field_value="empty_list"
            )
        
        with SafeFileOperations.safe_open(file_path, 'w') as handle:
            SeqIO.write(sequences, handle, "fasta")
        
        self.logger.info(f"Wrote {len(sequences)} sequences to {file_path}")
    
    def process(self, file_path: PathLike) -> Dict[str, Any]:
        """
        Process FASTA file and return metadata.
        
        Args:
            file_path: Path to FASTA file
            
        Returns:
            Dictionary with file statistics
        """
        sequences = self.read_file(file_path)
        
        total_length = sum(len(seq.seq) for seq in sequences)
        avg_length = total_length / len(sequences) if sequences else 0
        
        return {
            "file_path": str(file_path),
            "sequence_count": len(sequences),
            "total_length": total_length,
            "average_length": avg_length,
            "sequence_ids": [seq.id for seq in sequences],
        }
    
    def filter_sequences(
        self,
        file_path: PathLike,
        output_path: PathLike,
        min_length: Optional[int] = None,
        max_length: Optional[int] = None,
        ids_to_keep: Optional[List[str]] = None,
        ids_to_exclude: Optional[List[str]] = None,
    ) -> int:
        """
        Filter sequences based on criteria.
        
        Args:
            file_path: Input FASTA file
            output_path: Output FASTA file
            min_length: Minimum sequence length
            max_length: Maximum sequence length
            ids_to_keep: List of sequence IDs to keep
            ids_to_exclude: List of sequence IDs to exclude
            
        Returns:
            Number of sequences kept
        """
        sequences = self.read_file(file_path)
        filtered_sequences = []
        
        for seq in sequences:
            # Length filtering
            if min_length and len(seq.seq) < min_length:
                continue
            if max_length and len(seq.seq) > max_length:
                continue
            
            # ID filtering
            if ids_to_keep and seq.id not in ids_to_keep:
                continue
            if ids_to_exclude and seq.id in ids_to_exclude:
                continue
            
            filtered_sequences.append(seq)
        
        self.write_file(filtered_sequences, output_path)
        
        self.logger.info(
            f"Filtered {len(sequences)} -> {len(filtered_sequences)} sequences"
        )
        
        return len(filtered_sequences)


class GffProcessor(FileProcessor):
    """Processor for GFF/GTF format files."""
    
    def __init__(self, **kwargs):
        super().__init__(
            supported_formats=['.gff', '.gff3', '.gtf'],
            **kwargs
        )
    
    def validate_input(self, file_path: PathLike) -> None:
        """Validate GFF file."""
        super().validate_input(file_path)
        validate_gff_file(file_path)
    
    def read_file(self, file_path: PathLike) -> List[GeneInfo]:
        """
        Read GFF file and return gene information.
        
        Args:
            file_path: Path to GFF file
            
        Returns:
            List of gene information objects
        """
        self.validate_input(file_path)
        
        genes = []
        with SafeFileOperations.safe_open(file_path, 'r') as handle:
            for line_num, line in enumerate(handle, 1):
                line = line.strip()
                
                # Skip comments and empty lines
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) != 9:
                    continue
                
                # Parse GFF fields
                seqname, source, feature, start, end, score, strand, frame, attributes = parts
                
                # Only process gene features
                if feature.lower() != 'gene':
                    continue
                
                # Parse attributes
                attr_dict = self._parse_attributes(attributes)
                
                gene_info = GeneInfo(
                    gene_id=attr_dict.get('ID', f"gene_{line_num}"),
                    chromosome=seqname,
                    start=int(start),
                    end=int(end),
                    strand=strand,
                    gene_name=attr_dict.get('Name'),
                    description=attr_dict.get('Note'),
                    functional_annotation=attr_dict.get('function'),
                )
                
                genes.append(gene_info)
        
        self.logger.info(f"Read {len(genes)} genes from {file_path}")
        return genes
    
    def write_file(self, genes: List[GeneInfo], file_path: PathLike) -> None:
        """
        Write gene information to GFF file.
        
        Args:
            genes: List of gene information objects
            file_path: Output file path
        """
        with SafeFileOperations.safe_open(file_path, 'w') as handle:
            handle.write("##gff-version 3\n")
            
            for gene in genes:
                attributes = f"ID={gene.gene_id}"
                if gene.gene_name:
                    attributes += f";Name={gene.gene_name}"
                if gene.description:
                    attributes += f";Note={gene.description}"
                
                line = '\t'.join([
                    gene.chromosome,
                    'plant_mgc_analysis',
                    'gene',
                    str(gene.start),
                    str(gene.end),
                    '.',
                    gene.strand,
                    '.',
                    attributes
                ])
                
                handle.write(line + '\n')
        
        self.logger.info(f"Wrote {len(genes)} genes to {file_path}")
    
    def process(self, file_path: PathLike) -> Dict[str, Any]:
        """
        Process GFF file and return statistics.
        
        Args:
            file_path: Path to GFF file
            
        Returns:
            Dictionary with file statistics
        """
        genes = self.read_file(file_path)
        
        chromosomes = list(set(gene.chromosome for gene in genes))
        strand_counts = {'+': 0, '-': 0, '.': 0}
        
        for gene in genes:
            strand_counts[gene.strand] = strand_counts.get(gene.strand, 0) + 1
        
        return {
            "file_path": str(file_path),
            "gene_count": len(genes),
            "chromosome_count": len(chromosomes),
            "chromosomes": chromosomes,
            "strand_distribution": strand_counts,
        }
    
    def _parse_attributes(self, attributes: str) -> Dict[str, str]:
        """Parse GFF attributes field."""
        attr_dict = {}
        
        for attr in attributes.split(';'):
            if '=' in attr:
                key, value = attr.split('=', 1)
                attr_dict[key.strip()] = value.strip()
        
        return attr_dict


class CsvProcessor(FileProcessor):
    """Processor for CSV format files."""
    
    def __init__(self, **kwargs):
        super().__init__(
            supported_formats=['.csv', '.tsv'],
            **kwargs
        )
    
    def read_file(self, file_path: PathLike) -> pd.DataFrame:
        """
        Read CSV file and return DataFrame.
        
        Args:
            file_path: Path to CSV file
            
        Returns:
            Pandas DataFrame
        """
        self.validate_input(file_path)
        
        path = Path(file_path)
        separator = '\t' if path.suffix == '.tsv' else ','
        
        try:
            with SafeFileOperations.safe_open(file_path, 'r') as handle:
                df = pd.read_csv(handle, sep=separator)
            
            self.logger.info(f"Read CSV with shape {df.shape} from {file_path}")
            return df
            
        except Exception as e:
            raise FileSystemError(
                f"Error reading CSV file: {e}",
                file_path=str(file_path),
                operation="read_csv"
            ) from e
    
    def write_file(self, data: pd.DataFrame, file_path: PathLike) -> None:
        """
        Write DataFrame to CSV file.
        
        Args:
            data: DataFrame to write
            file_path: Output file path
        """
        path = Path(file_path)
        separator = '\t' if path.suffix == '.tsv' else ','
        
        with SafeFileOperations.safe_open(file_path, 'w') as handle:
            data.to_csv(handle, sep=separator, index=False)
        
        self.logger.info(f"Wrote CSV with shape {data.shape} to {file_path}")
    
    def process(self, file_path: PathLike) -> Dict[str, Any]:
        """
        Process CSV file and return metadata.
        
        Args:
            file_path: Path to CSV file
            
        Returns:
            Dictionary with file statistics
        """
        df = self.read_file(file_path)
        
        return {
            "file_path": str(file_path),
            "row_count": len(df),
            "column_count": len(df.columns),
            "columns": list(df.columns),
            "dtypes": df.dtypes.to_dict(),
            "null_counts": df.isnull().sum().to_dict(),
        }


class BlastResultProcessor(FileProcessor):
    """Processor for BLAST result files."""
    
    def __init__(self, **kwargs):
        super().__init__(
            supported_formats=['.xml', '.txt', '.json'],
            **kwargs
        )
    
    def read_file(self, file_path: PathLike) -> List[Dict[str, Any]]:
        """
        Read BLAST results file.
        
        Args:
            file_path: Path to BLAST results file
            
        Returns:
            List of BLAST hit dictionaries
        """
        self.validate_input(file_path)
        
        path = Path(file_path)
        results = []
        
        try:
            if path.suffix == '.xml':
                results = self._read_blast_xml(path)
            elif path.suffix == '.txt':
                results = self._read_blast_tabular(path)
            elif path.suffix == '.json':
                results = self._read_blast_json(path)
            else:
                raise ValidationError(
                    f"Unsupported BLAST format: {path.suffix}",
                    field_name="file_format",
                    field_value=path.suffix
                )
            
            self.logger.info(f"Read {len(results)} BLAST results from {file_path}")
            return results
            
        except Exception as e:
            raise FileSystemError(
                f"Error reading BLAST file: {e}",
                file_path=str(file_path),
                operation="read_blast"
            ) from e
    
    def write_file(self, results: List[Dict[str, Any]], file_path: PathLike) -> None:
        """
        Write BLAST results to file.
        
        Args:
            results: List of BLAST result dictionaries
            file_path: Output file path
        """
        path = Path(file_path)
        
        if path.suffix == '.json':
            with SafeFileOperations.safe_open(file_path, 'w') as handle:
                json.dump(results, handle, indent=2, default=str)
        elif path.suffix == '.csv':
            df = pd.DataFrame(results)
            CsvProcessor().write_file(df, file_path)
        else:
            raise ValidationError(
                f"Unsupported output format: {path.suffix}",
                field_name="output_format",
                field_value=path.suffix
            )
        
        self.logger.info(f"Wrote {len(results)} BLAST results to {file_path}")
    
    def process(self, file_path: PathLike) -> Dict[str, Any]:
        """
        Process BLAST file and return statistics.
        
        Args:
            file_path: Path to BLAST file
            
        Returns:
            Dictionary with statistics
        """
        results = self.read_file(file_path)
        
        if not results:
            return {"file_path": str(file_path), "hit_count": 0}
        
        evalues = [hit.get('evalue', float('inf')) for hit in results]
        bit_scores = [hit.get('bit_score', 0) for hit in results]
        identities = [hit.get('identity', 0) for hit in results]
        
        return {
            "file_path": str(file_path),
            "hit_count": len(results),
            "mean_evalue": sum(evalues) / len(evalues),
            "mean_bit_score": sum(bit_scores) / len(bit_scores),
            "mean_identity": sum(identities) / len(identities),
            "unique_queries": len(set(hit.get('query_id') for hit in results)),
            "unique_subjects": len(set(hit.get('subject_id') for hit in results)),
        }
    
    def _read_blast_xml(self, file_path: Path) -> List[Dict[str, Any]]:
        """Read BLAST XML format."""
        results = []
        
        with SafeFileOperations.safe_open(file_path, 'r') as handle:
            for blast_record in SearchIO.parse(handle, 'blast-xml'):
                for hit in blast_record:
                    for hsp in hit:
                        result = {
                            'query_id': blast_record.id,
                            'subject_id': hit.id,
                            'identity': hsp.ident_num / hsp.aln_len * 100,
                            'alignment_length': hsp.aln_len,
                            'mismatches': hsp.aln_len - hsp.ident_num,
                            'gap_opens': hsp.gap_num,
                            'query_start': hsp.query_start,
                            'query_end': hsp.query_end,
                            'subject_start': hsp.hit_start,
                            'subject_end': hsp.hit_end,
                            'evalue': hsp.evalue,
                            'bit_score': hsp.bitscore,
                        }
                        results.append(result)
        
        return results
    
    def _read_blast_tabular(self, file_path: Path) -> List[Dict[str, Any]]:
        """Read BLAST tabular format."""
        results = []
        
        # Standard BLAST tabular columns
        columns = [
            'query_id', 'subject_id', 'identity', 'alignment_length',
            'mismatches', 'gap_opens', 'query_start', 'query_end',
            'subject_start', 'subject_end', 'evalue', 'bit_score'
        ]
        
        with SafeFileOperations.safe_open(file_path, 'r') as handle:
            reader = csv.DictReader(handle, fieldnames=columns, delimiter='\t')
            
            for row in reader:
                # Skip comment lines
                if row['query_id'].startswith('#'):
                    continue
                
                # Convert numeric fields
                for field in ['identity', 'alignment_length', 'mismatches', 
                             'gap_opens', 'query_start', 'query_end',
                             'subject_start', 'subject_end', 'evalue', 'bit_score']:
                    try:
                        if field in ['evalue', 'bit_score', 'identity']:
                            row[field] = float(row[field])
                        else:
                            row[field] = int(row[field])
                    except (ValueError, TypeError):
                        row[field] = None
                
                results.append(row)
        
        return results
    
    def _read_blast_json(self, file_path: Path) -> List[Dict[str, Any]]:
        """Read BLAST JSON format."""
        with SafeFileOperations.safe_open(file_path, 'r') as handle:
            data = json.load(handle)
        
        if isinstance(data, list):
            return data
        elif isinstance(data, dict) and 'results' in data:
            return data['results']
        else:
            raise ValidationError(
                "Invalid BLAST JSON format",
                field_name="json_structure",
                field_value=type(data).__name__
            )


def read_fasta(file_path: PathLike) -> List[SeqRecord]:
    """Convenience function to read FASTA file."""
    processor = FastaProcessor()
    return processor.read_file(file_path)


def write_fasta(sequences: List[SeqRecord], file_path: PathLike) -> None:
    """Convenience function to write FASTA file."""
    processor = FastaProcessor()
    processor.write_file(sequences, file_path)


def read_gff(file_path: PathLike) -> List[GeneInfo]:
    """Convenience function to read GFF file."""
    processor = GffProcessor()
    return processor.read_file(file_path)


def safe_file_operation(operation: callable, *args, **kwargs) -> Any:
    """
    Safely execute file operation with error handling.
    
    Args:
        operation: File operation function
        *args: Operation arguments
        **kwargs: Operation keyword arguments
        
    Returns:
        Operation result
        
    Raises:
        FileSystemError: If operation fails
    """
    try:
        return operation(*args, **kwargs)
    except Exception as e:
        raise FileSystemError(
            f"File operation failed: {e}",
            operation=operation.__name__
        ) from e