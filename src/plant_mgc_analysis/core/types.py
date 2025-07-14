"""
Type definitions for Plant MGC Analysis Pipeline.

This module defines the core data structures and types used throughout
the application for better type safety and documentation.
"""

from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass, field
from pathlib import Path
from datetime import datetime
from enum import Enum


class AnalysisType(str, Enum):
    """Types of analyses supported by the pipeline."""
    
    SLIDING_WINDOW = "sliding_window"
    BLAST_SEARCH = "blast_search"
    PHYLOGENETIC = "phylogenetic"
    MACHINE_LEARNING = "machine_learning"
    PATHWAY_ANALYSIS = "pathway_analysis"
    COMPARATIVE_GENOMICS = "comparative_genomics"


class MGCStatus(str, Enum):
    """Status of MGC candidate."""
    
    CANDIDATE = "candidate"
    CONFIRMED = "confirmed"
    REJECTED = "rejected"
    UNKNOWN = "unknown"


class DatabaseType(str, Enum):
    """Supported database types."""
    
    KEGG = "kegg"
    MIBIG = "mibig"
    PLANTCYC = "plantcyc"
    ENSEMBL = "ensembl"
    PHYTOZOME = "phytozome"
    PLAZA = "plaza"


@dataclass
class GenomeData:
    """Represents genome data and metadata."""
    
    organism_name: str
    genome_file: Path
    annotation_file: Optional[Path] = None
    protein_file: Optional[Path] = None
    species: Optional[str] = None
    strain: Optional[str] = None
    assembly_version: Optional[str] = None
    database_id: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Validate genome data after initialization."""
        if not self.genome_file.exists():
            raise FileNotFoundError(f"Genome file not found: {self.genome_file}")
        
        if self.annotation_file and not self.annotation_file.exists():
            raise FileNotFoundError(f"Annotation file not found: {self.annotation_file}")


@dataclass
class GeneInfo:
    """Represents gene information."""
    
    gene_id: str
    chromosome: str
    start: int
    end: int
    strand: str
    gene_name: Optional[str] = None
    protein_id: Optional[str] = None
    description: Optional[str] = None
    functional_annotation: Optional[str] = None
    pathway_ids: List[str] = field(default_factory=list)
    metabolic_role: Optional[str] = None
    
    @property
    def length(self) -> int:
        """Return gene length."""
        return abs(self.end - self.start)
    
    @property
    def location(self) -> str:
        """Return formatted location string."""
        return f"{self.chromosome}:{self.start}-{self.end}({self.strand})"


@dataclass
class MGCCandidate:
    """Represents a metabolic gene cluster candidate."""
    
    cluster_id: str
    chromosome: str
    start: int
    end: int
    genes: List[GeneInfo]
    organism: str
    analysis_method: str
    confidence_score: float
    status: MGCStatus = MGCStatus.CANDIDATE
    pathway_id: Optional[str] = None
    pathway_name: Optional[str] = None
    metabolic_class: Optional[str] = None
    validation_data: Dict[str, Any] = field(default_factory=dict)
    created_at: datetime = field(default_factory=datetime.now)
    
    @property
    def length(self) -> int:
        """Return cluster length in base pairs."""
        return abs(self.end - self.start)
    
    @property
    def gene_count(self) -> int:
        """Return number of genes in cluster."""
        return len(self.genes)
    
    @property
    def location(self) -> str:
        """Return formatted location string."""
        return f"{self.chromosome}:{self.start}-{self.end}"
    
    @property
    def gene_density(self) -> float:
        """Return gene density (genes per kb)."""
        if self.length == 0:
            return 0.0
        return (self.gene_count / self.length) * 1000
    
    def get_metabolic_genes(self) -> List[GeneInfo]:
        """Return genes with metabolic annotations."""
        return [gene for gene in self.genes if gene.metabolic_role is not None]
    
    def get_pathway_genes(self, pathway_id: str) -> List[GeneInfo]:
        """Return genes associated with a specific pathway."""
        return [gene for gene in self.genes if pathway_id in gene.pathway_ids]


@dataclass
class BlastResult:
    """Represents BLAST search result."""
    
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
    query_coverage: Optional[float] = None
    subject_coverage: Optional[float] = None
    
    @property
    def query_length(self) -> int:
        """Return query alignment length."""
        return abs(self.query_end - self.query_start)
    
    @property
    def subject_length(self) -> int:
        """Return subject alignment length."""
        return abs(self.subject_end - self.subject_start)


@dataclass
class AnalysisResult:
    """Represents analysis results."""
    
    analysis_id: str
    analysis_type: AnalysisType
    input_data: Dict[str, Any]
    results: Dict[str, Any]
    parameters: Dict[str, Any]
    status: str
    error_message: Optional[str] = None
    execution_time: Optional[float] = None
    created_at: datetime = field(default_factory=datetime.now)
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def success(self) -> bool:
        """Return True if analysis was successful."""
        return self.status == "success"
    
    def get_mgc_candidates(self) -> List[MGCCandidate]:
        """Extract MGC candidates from results."""
        if "mgc_candidates" in self.results:
            return self.results["mgc_candidates"]
        return []
    
    def get_statistics(self) -> Dict[str, Any]:
        """Extract statistical results."""
        return self.results.get("statistics", {})


@dataclass
class PhylogeneticData:
    """Represents phylogenetic analysis data."""
    
    tree_file: Path
    alignment_file: Path
    sequence_ids: List[str]
    tree_format: str = "newick"
    alignment_format: str = "fasta"
    bootstrap_support: Optional[List[float]] = None
    substitution_model: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PipelineConfig:
    """Configuration for analysis pipeline."""
    
    pipeline_name: str
    analysis_types: List[AnalysisType]
    input_files: List[Path]
    output_dir: Path
    parameters: Dict[str, Any] = field(default_factory=dict)
    parallel_jobs: int = 1
    use_cache: bool = True
    cache_dir: Optional[Path] = None
    
    def validate(self) -> None:
        """Validate pipeline configuration."""
        if not self.input_files:
            raise ValueError("No input files specified")
        
        for file_path in self.input_files:
            if not file_path.exists():
                raise FileNotFoundError(f"Input file not found: {file_path}")
        
        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True, exist_ok=True)


# Type aliases for commonly used types
PathLike = Union[str, Path]
BlastResults = List[BlastResult]
MGCCandidates = List[MGCCandidate]
AnalysisResults = List[AnalysisResult]
GenomeDataList = List[GenomeData]
GeneInfoList = List[GeneInfo]

# Database query types
DatabaseQuery = Dict[str, Any]
DatabaseResult = Dict[str, Any]
DatabaseResults = List[DatabaseResult]

# Statistical analysis types
StatisticalTest = str
PValue = float
StatisticalResults = Dict[str, Union[PValue, str, float]]

# Machine learning types
MLModel = Any  # Will be more specific when ML module is implemented
MLFeatures = Dict[str, Any]
MLPrediction = Dict[str, Any]
MLResults = List[MLPrediction]