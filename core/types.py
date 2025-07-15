"""
Type definitions and data structures for bioinformatics operations.

This module provides common data types, enums, and data classes used throughout
the bioinformatics toolkit.
"""

from typing import Any, Dict, List, Optional, Union, Tuple
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
import datetime


class ProcessingStatus(Enum):
    """Status of processing operations."""
    PENDING = "pending"
    IN_PROGRESS = "in_progress"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class FileFormat(Enum):
    """Supported file formats."""
    FASTA = "fasta"
    CSV = "csv"
    TSV = "tsv"
    JSON = "json"
    XML = "xml"
    GFF = "gff"
    GFF3 = "gff3"
    GTF = "gtf"
    GENBANK = "genbank"
    BLAST_XML = "blast_xml"
    BLAST_TSV = "blast_tsv"


class SequenceType(Enum):
    """Types of biological sequences."""
    DNA = "dna"
    RNA = "rna"
    PROTEIN = "protein"
    NUCLEOTIDE = "nucleotide"


class DatabaseType(Enum):
    """Supported database types."""
    KEGG = "kegg"
    NCBI = "ncbi"
    UNIPROT = "uniprot"
    MIBIG = "mibig"
    PLANTCYC = "plantcyc"
    ENSEMBL = "ensembl"
    PHYTOZOME = "phytozome"


class AnalysisType(Enum):
    """Types of bioinformatics analyses."""
    BLAST_SEARCH = "blast_search"
    SLIDING_WINDOW = "sliding_window"
    SEQUENCE_ALIGNMENT = "sequence_alignment"
    PHYLOGENETIC = "phylogenetic"
    STATISTICAL = "statistical"
    MGC_DETECTION = "mgc_detection"


@dataclass
class Coordinates:
    """Genomic coordinates."""
    chromosome: str
    start: int
    end: int
    strand: Optional[str] = None
    
    def __post_init__(self):
        if self.start > self.end:
            raise ValueError(f"Start position {self.start} cannot be greater than end position {self.end}")
    
    @property
    def length(self) -> int:
        """Length of the genomic region."""
        return self.end - self.start + 1


@dataclass
class GeneInfo:
    """Information about a gene."""
    gene_id: str
    coordinates: Coordinates
    gene_name: Optional[str] = None
    function: Optional[str] = None
    ec_number: Optional[str] = None
    pathway: Optional[str] = None
    sequence: Optional[str] = None
    protein_sequence: Optional[str] = None


@dataclass
class SequenceData:
    """Biological sequence data."""
    sequence_id: str
    sequence: str
    sequence_type: SequenceType
    description: Optional[str] = None
    organism: Optional[str] = None
    length: Optional[int] = field(init=False)
    
    def __post_init__(self):
        self.length = len(self.sequence)


@dataclass
class GenomeData:
    """Genome information and data."""
    organism_name: str
    genome_file: Path
    annotation_file: Optional[Path] = None
    genes: List[GeneInfo] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def gene_count(self) -> int:
        """Number of genes in the genome."""
        return len(self.genes)


@dataclass
class BlastHit:
    """BLAST search result hit."""
    query_id: str
    subject_id: str
    identity: float
    alignment_length: int
    mismatches: int
    gap_opens: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    evalue: float
    bit_score: float
    subject_description: Optional[str] = None


@dataclass
class AnalysisResult:
    """Result of a bioinformatics analysis."""
    analysis_type: AnalysisType
    status: ProcessingStatus
    results: Dict[str, Any]
    metadata: Dict[str, Any] = field(default_factory=dict)
    timestamp: datetime.datetime = field(default_factory=datetime.datetime.now)
    processing_time: Optional[float] = None
    error_message: Optional[str] = None
    
    def is_successful(self) -> bool:
        """Check if analysis completed successfully."""
        return self.status == ProcessingStatus.COMPLETED


@dataclass
class MGCCandidate:
    """Metabolic gene cluster candidate."""
    cluster_id: str
    coordinates: Coordinates
    genes: List[GeneInfo]
    confidence_score: float
    pathway_annotation: Optional[str] = None
    ec_numbers: List[str] = field(default_factory=list)
    status: ProcessingStatus = ProcessingStatus.PENDING
    
    @property
    def gene_count(self) -> int:
        """Number of genes in the cluster."""
        return len(self.genes)


@dataclass
class ProcessingConfig:
    """Configuration for processing operations."""
    batch_size: int = 100
    max_workers: int = 4
    timeout: int = 300
    retry_attempts: int = 3
    use_cache: bool = True
    output_format: FileFormat = FileFormat.CSV
    

@dataclass
class DatabaseQuery:
    """Database query parameters."""
    database: DatabaseType
    query: str
    limit: Optional[int] = None
    filters: Dict[str, Any] = field(default_factory=dict)
    use_cache: bool = True


# Type aliases for common combinations
FilePath = Union[str, Path]
Numeric = Union[int, float]
QueryResults = List[Dict[str, Any]]
ParameterDict = Dict[str, Any]