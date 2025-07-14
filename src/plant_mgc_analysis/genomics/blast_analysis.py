"""
BLAST Analysis for Plant MGC Analysis Pipeline.

This module provides object-oriented BLAST analysis functionality with proper
database management, result processing, and statistical validation.
"""

import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple, Iterator, Any
from dataclasses import dataclass, field
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import json

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Blast import NCBIXML, NCBIWWW
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from ..core.base import AnalysisEngine, DatabaseInterface
from ..core.types import (
    GenomeData,
    AnalysisResult,
    AnalysisType,
    BlastResult,
    DatabaseType,
)
from ..core.exceptions import (
    AnalysisError,
    DatabaseError,
    ValidationError,
    ComputeError,
)
from ..utils.file_operations import FastaProcessor, CsvProcessor, SafeFileOperations
from ..utils.validation import validate_analysis_parameters


@dataclass
class BlastSearchParameters:
    """Parameters for BLAST search."""
    
    program: str = "blastp"  # blastp, blastn, blastx, tblastn, tblastx
    evalue: float = 1e-5
    max_target_seqs: int = 100
    num_threads: int = 4
    word_size: Optional[int] = None
    matrix: str = "BLOSUM62"
    gap_costs: Tuple[int, int] = (11, 1)  # gap open, gap extend
    seg_filter: bool = True
    dust_filter: bool = True
    
    def to_blast_args(self) -> List[str]:
        """Convert parameters to BLAST command line arguments."""
        args = [
            "-evalue", str(self.evalue),
            "-max_target_seqs", str(self.max_target_seqs),
            "-num_threads", str(self.num_threads),
            "-matrix", self.matrix,
            "-gapopen", str(self.gap_costs[0]),
            "-gapextend", str(self.gap_costs[1]),
        ]
        
        if self.word_size:
            args.extend(["-word_size", str(self.word_size)])
        
        if not self.seg_filter:
            args.append("-seg no")
        
        if not self.dust_filter:
            args.append("-dust no")
        
        return args


@dataclass 
class BlastHit:
    """Represents a BLAST search hit."""
    
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
    
    # Additional computed fields
    query_coverage: Optional[float] = None
    subject_coverage: Optional[float] = None
    combined_score: Optional[float] = None
    
    # Source information
    source_file: Optional[str] = None
    subject_index: Optional[int] = None
    
    def calculate_coverage(self, query_length: int, subject_length: int) -> None:
        """Calculate query and subject coverage."""
        self.query_coverage = ((self.query_end - self.query_start + 1) / query_length) * 100
        self.subject_coverage = ((self.subject_end - self.subject_start + 1) / subject_length) * 100
    
    def calculate_combined_score(self, weights: Dict[str, float] = None) -> float:
        """
        Calculate combined score from multiple metrics.
        
        Args:
            weights: Weights for different metrics (bit_score, identity, coverage)
            
        Returns:
            Combined score
        """
        if weights is None:
            weights = {"bit_score": 0.4, "identity": 0.3, "coverage": 0.3}
        
        coverage = self.query_coverage or 0.0
        
        self.combined_score = (
            (self.bit_score * weights.get("bit_score", 0.4)) +
            (self.identity * weights.get("identity", 0.3)) +
            (coverage * weights.get("coverage", 0.3))
        )
        
        return self.combined_score
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "query_id": self.query_id,
            "subject_id": self.subject_id,
            "identity": self.identity,
            "alignment_length": self.alignment_length,
            "mismatches": self.mismatches,
            "gap_opens": self.gap_opens,
            "query_start": self.query_start,
            "query_end": self.query_end,
            "subject_start": self.subject_start,
            "subject_end": self.subject_end,
            "evalue": self.evalue,
            "bit_score": self.bit_score,
            "query_coverage": self.query_coverage,
            "subject_coverage": self.subject_coverage,
            "combined_score": self.combined_score,
            "source_file": self.source_file,
            "subject_index": self.subject_index,
        }


class BlastDatabase:
    """Manages BLAST database creation and operations."""
    
    def __init__(self, db_path: Path, db_type: str = "prot"):
        """
        Initialize BLAST database.
        
        Args:
            db_path: Path to database files (without extension)
            db_type: Database type ('nucl' or 'prot')
        """
        self.db_path = Path(db_path)
        self.db_type = db_type
        self._is_built = False
    
    def build_from_fasta(self, fasta_file: Path, title: Optional[str] = None) -> None:
        """
        Build BLAST database from FASTA file.
        
        Args:
            fasta_file: Input FASTA file
            title: Optional database title
            
        Raises:
            ComputeError: If database building fails
        """
        if self.exists():
            logger.info(f"BLAST database already exists: {self.db_path}")
            self._is_built = True
            return
        
        # Ensure output directory exists
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        
        cmd = [
            "makeblastdb",
            "-in", str(fasta_file),
            "-dbtype", self.db_type,
            "-out", str(self.db_path),
        ]
        
        if title:
            cmd.extend(["-title", title])
        
        try:
            logger.info(f"Building BLAST database: {self.db_path}")
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            logger.info("BLAST database built successfully")
            self._is_built = True
            
        except subprocess.CalledProcessError as e:
            raise ComputeError(
                f"Failed to build BLAST database: {e.stderr}",
                computation_type="makeblastdb",
                resource_type="database"
            ) from e
        except FileNotFoundError:
            raise ComputeError(
                "makeblastdb command not found. Please install BLAST+",
                computation_type="makeblastdb",
                resource_type="software"
            )
    
    def build_from_sequences(
        self,
        sequences: List[SeqRecord],
        title: Optional[str] = None
    ) -> None:
        """
        Build BLAST database from sequence records.
        
        Args:
            sequences: List of sequence records
            title: Optional database title
        """
        # Create temporary FASTA file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_file:
            SeqIO.write(sequences, tmp_file, "fasta")
            tmp_fasta = Path(tmp_file.name)
        
        try:
            self.build_from_fasta(tmp_fasta, title)
        finally:
            tmp_fasta.unlink()  # Clean up temporary file
    
    def exists(self) -> bool:
        """Check if database files exist."""
        extensions = ['.pin', '.phr', '.psq'] if self.db_type == 'prot' else ['.nin', '.nhr', '.nsq']
        
        # Check for single-volume database
        if all((self.db_path.parent / f"{self.db_path.name}{ext}").exists() for ext in extensions):
            return True
        
        # Check for multi-volume database
        if any((self.db_path.parent / f"{self.db_path.name}.00{ext}").exists() for ext in extensions):
            return True
        
        return False
    
    def get_info(self) -> Dict[str, Any]:
        """Get database information using blastdbcmd."""
        if not self.exists():
            return {}
        
        try:
            result = subprocess.run(
                ["blastdbcmd", "-db", str(self.db_path), "-info"],
                capture_output=True,
                text=True,
                check=True
            )
            
            # Parse output (simplified)
            info = {"raw_output": result.stdout}
            return info
            
        except subprocess.CalledProcessError:
            return {}


class BlastSearchEngine:
    """Executes BLAST searches with various configurations."""
    
    def __init__(self, database: BlastDatabase):
        """
        Initialize BLAST search engine.
        
        Args:
            database: BLAST database to search against
        """
        self.database = database
        self.temp_files: List[Path] = []
    
    def search(
        self,
        query_sequences: List[SeqRecord],
        parameters: BlastSearchParameters,
        output_format: str = "xml"
    ) -> Path:
        """
        Perform BLAST search.
        
        Args:
            query_sequences: Query sequences
            parameters: BLAST parameters
            output_format: Output format (xml, json, tabular)
            
        Returns:
            Path to results file
            
        Raises:
            ComputeError: If BLAST search fails
        """
        if not self.database.exists():
            raise DatabaseError(
                "BLAST database does not exist",
                database_name=str(self.database.db_path),
                operation="blast_search"
            )
        
        # Create temporary query file
        query_file = Path(tempfile.mktemp(suffix='.fasta'))
        with SafeFileOperations.safe_open(query_file, 'w') as handle:
            SeqIO.write(query_sequences, handle, "fasta")
        self.temp_files.append(query_file)
        
        # Create output file
        output_file = Path(tempfile.mktemp(suffix=f'.{output_format}'))
        self.temp_files.append(output_file)
        
        # Build command
        cmd = [
            parameters.program,
            "-query", str(query_file),
            "-db", str(self.database.db_path),
            "-out", str(output_file),
        ]
        
        # Add output format
        if output_format == "xml":
            cmd.extend(["-outfmt", "5"])
        elif output_format == "json":
            cmd.extend(["-outfmt", "15"])
        elif output_format == "tabular":
            cmd.extend(["-outfmt", "6"])
        
        # Add other parameters
        cmd.extend(parameters.to_blast_args())
        
        try:
            logger.info(f"Running BLAST search with {len(query_sequences)} queries")
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            logger.info(f"BLAST search completed: {output_file}")
            return output_file
            
        except subprocess.CalledProcessError as e:
            raise ComputeError(
                f"BLAST search failed: {e.stderr}",
                computation_type="blast_search",
                resource_type="computation"
            ) from e
    
    def cleanup(self) -> None:
        """Clean up temporary files."""
        for temp_file in self.temp_files:
            if temp_file.exists():
                temp_file.unlink()
        self.temp_files.clear()


class BlastResultParser:
    """Parses BLAST results from various formats."""
    
    @staticmethod
    def parse_xml(xml_file: Path) -> List[BlastHit]:
        """
        Parse BLAST XML results.
        
        Args:
            xml_file: Path to XML results file
            
        Returns:
            List of BLAST hits
        """
        hits = []
        
        with SafeFileOperations.safe_open(xml_file, 'r') as handle:
            blast_records = NCBIXML.parse(handle)
            
            for record in blast_records:
                query_id = record.query
                query_length = record.query_length
                
                for alignment in record.alignments:
                    subject_id = alignment.hit_def
                    subject_length = alignment.length
                    
                    for hsp in alignment.hsps:
                        hit = BlastHit(
                            query_id=query_id,
                            subject_id=subject_id,
                            identity=(hsp.identities / hsp.align_length) * 100,
                            alignment_length=hsp.align_length,
                            mismatches=hsp.align_length - hsp.identities,
                            gap_opens=hsp.gaps,
                            query_start=hsp.query_start,
                            query_end=hsp.query_end,
                            subject_start=hsp.sbjct_start,
                            subject_end=hsp.sbjct_end,
                            evalue=hsp.expect,
                            bit_score=hsp.bits,
                        )
                        
                        hit.calculate_coverage(query_length, subject_length)
                        hits.append(hit)
        
        return hits
    
    @staticmethod
    def parse_tabular(tabular_file: Path) -> List[BlastHit]:
        """
        Parse BLAST tabular results.
        
        Args:
            tabular_file: Path to tabular results file
            
        Returns:
            List of BLAST hits
        """
        hits = []
        
        columns = [
            'query_id', 'subject_id', 'identity', 'alignment_length',
            'mismatches', 'gap_opens', 'query_start', 'query_end',
            'subject_start', 'subject_end', 'evalue', 'bit_score'
        ]
        
        processor = CsvProcessor()
        df = pd.read_csv(tabular_file, sep='\t', names=columns, comment='#')
        
        for _, row in df.iterrows():
            hit = BlastHit(
                query_id=row['query_id'],
                subject_id=row['subject_id'],
                identity=float(row['identity']),
                alignment_length=int(row['alignment_length']),
                mismatches=int(row['mismatches']),
                gap_opens=int(row['gap_opens']),
                query_start=int(row['query_start']),
                query_end=int(row['query_end']),
                subject_start=int(row['subject_start']),
                subject_end=int(row['subject_end']),
                evalue=float(row['evalue']),
                bit_score=float(row['bit_score']),
            )
            hits.append(hit)
        
        return hits
    
    @staticmethod
    def parse_json(json_file: Path) -> List[BlastHit]:
        """
        Parse BLAST JSON results.
        
        Args:
            json_file: Path to JSON results file
            
        Returns:
            List of BLAST hits
        """
        with SafeFileOperations.safe_open(json_file, 'r') as handle:
            data = json.load(handle)
        
        hits = []
        
        # Parse JSON structure (simplified implementation)
        if 'BlastOutput2' in data:
            for result in data['BlastOutput2']:
                # Implementation would depend on BLAST JSON format
                # This is a placeholder for the actual parsing logic
                pass
        
        return hits


class BlastAnalyzer(AnalysisEngine):
    """
    BLAST-based analysis engine for metabolic gene cluster detection.
    
    This class provides comprehensive BLAST analysis functionality including
    database management, search execution, and result processing.
    """
    
    def __init__(self, **kwargs):
        """Initialize BLAST analyzer."""
        super().__init__(
            analysis_type=AnalysisType.BLAST_SEARCH,
            **kwargs
        )
        
        self.database: Optional[BlastDatabase] = None
        self.search_engine: Optional[BlastSearchEngine] = None
        self.reference_databases: Dict[str, Path] = {}
    
    def validate_input(self, genome_data: GenomeData) -> None:
        """Validate input genome data."""
        if not genome_data.protein_file:
            raise ValidationError(
                "Protein file required for BLAST analysis",
                field_name="protein_file",
                field_value=None
            )
        
        if not genome_data.protein_file.exists():
            raise ValidationError(
                f"Protein file not found: {genome_data.protein_file}",
                field_name="protein_file",
                field_value=str(genome_data.protein_file)
            )
    
    def validate_parameters(self, parameters: Dict[str, Any]) -> None:
        """Validate BLAST analysis parameters."""
        super().validate_parameters(parameters)
        
        validated = validate_analysis_parameters("blast_search", parameters)
        
        # Validate BLAST-specific parameters
        if 'program' in parameters:
            valid_programs = ['blastp', 'blastn', 'blastx', 'tblastn', 'tblastx']
            if parameters['program'] not in valid_programs:
                raise ValidationError(
                    f"Invalid BLAST program: {parameters['program']}. "
                    f"Valid options: {valid_programs}",
                    field_name="program",
                    field_value=parameters['program']
                )
    
    def load_reference_database(self, database_name: str, database_path: Path) -> None:
        """
        Load reference database for BLAST searches.
        
        Args:
            database_name: Name of the database
            database_path: Path to database files
        """
        self.reference_databases[database_name] = database_path
        logger.info(f"Loaded reference database: {database_name} at {database_path}")
    
    def create_genome_database(self, genome_data: GenomeData) -> BlastDatabase:
        """
        Create BLAST database from genome data.
        
        Args:
            genome_data: Genome data
            
        Returns:
            BLAST database object
        """
        db_path = self.temp_dir / f"{genome_data.organism_name}_blastdb"
        database = BlastDatabase(db_path, "prot")
        
        # Load protein sequences
        processor = FastaProcessor()
        sequences = processor.read_file(genome_data.protein_file)
        
        # Build database
        database.build_from_sequences(sequences, f"{genome_data.organism_name} proteins")
        
        return database
    
    def search_against_reference(
        self,
        query_sequences: List[SeqRecord],
        reference_db: str,
        parameters: BlastSearchParameters,
    ) -> List[BlastHit]:
        """
        Search sequences against reference database.
        
        Args:
            query_sequences: Query sequences
            reference_db: Reference database name
            parameters: BLAST parameters
            
        Returns:
            List of BLAST hits
        """
        if reference_db not in self.reference_databases:
            raise ValidationError(
                f"Reference database not found: {reference_db}",
                field_name="reference_database",
                field_value=reference_db
            )
        
        # Create database object
        db_path = self.reference_databases[reference_db]
        database = BlastDatabase(db_path, "prot")
        
        # Create search engine
        search_engine = BlastSearchEngine(database)
        
        try:
            # Perform search
            results_file = search_engine.search(query_sequences, parameters, "xml")
            
            # Parse results
            hits = BlastResultParser.parse_xml(results_file)
            
            # Calculate combined scores
            for hit in hits:
                hit.calculate_combined_score()
            
            logger.info(f"Found {len(hits)} BLAST hits against {reference_db}")
            return hits
            
        finally:
            search_engine.cleanup()
    
    def filter_hits_by_quality(
        self,
        hits: List[BlastHit],
        min_identity: float = 30.0,
        min_coverage: float = 50.0,
        max_evalue: float = 1e-5,
    ) -> List[BlastHit]:
        """
        Filter BLAST hits by quality criteria.
        
        Args:
            hits: List of BLAST hits
            min_identity: Minimum identity percentage
            min_coverage: Minimum query coverage percentage
            max_evalue: Maximum E-value
            
        Returns:
            Filtered hits
        """
        filtered_hits = []
        
        for hit in hits:
            if (hit.identity >= min_identity and
                hit.evalue <= max_evalue and
                (hit.query_coverage or 0) >= min_coverage):
                filtered_hits.append(hit)
        
        logger.info(f"Filtered {len(hits)} -> {len(filtered_hits)} hits")
        return filtered_hits
    
    def cluster_hits_by_region(
        self,
        hits: List[BlastHit],
        max_distance: int = 100000,  # 100kb
    ) -> List[List[BlastHit]]:
        """
        Cluster hits by genomic region.
        
        Args:
            hits: List of BLAST hits
            max_distance: Maximum distance between hits in a cluster
            
        Returns:
            List of hit clusters
        """
        # Group hits by subject sequence
        subject_hits = {}
        for hit in hits:
            if hit.subject_id not in subject_hits:
                subject_hits[hit.subject_id] = []
            subject_hits[hit.subject_id].append(hit)
        
        clusters = []
        
        for subject_id, subject_hit_list in subject_hits.items():
            # Sort by position
            subject_hit_list.sort(key=lambda h: h.subject_start)
            
            current_cluster = [subject_hit_list[0]]
            
            for hit in subject_hit_list[1:]:
                # Check distance to last hit in current cluster
                last_hit = current_cluster[-1]
                distance = hit.subject_start - last_hit.subject_end
                
                if distance <= max_distance:
                    current_cluster.append(hit)
                else:
                    # Start new cluster
                    if len(current_cluster) > 1:  # Only keep multi-gene clusters
                        clusters.append(current_cluster)
                    current_cluster = [hit]
            
            # Add final cluster
            if len(current_cluster) > 1:
                clusters.append(current_cluster)
        
        logger.info(f"Found {len(clusters)} hit clusters")
        return clusters
    
    def analyze(
        self,
        genome_data: GenomeData,
        parameters: Optional[Dict[str, Any]] = None,
    ) -> AnalysisResult:
        """
        Perform BLAST analysis.
        
        Args:
            genome_data: Genome data to analyze
            parameters: Optional analysis parameters
            
        Returns:
            Analysis result
        """
        if parameters:
            self.validate_parameters(parameters)
        
        # Extract parameters
        reference_db = parameters.get("reference_database", "mibig")
        blast_params = BlastSearchParameters(
            evalue=parameters.get("evalue", 1e-5),
            max_target_seqs=parameters.get("max_targets", 100),
            num_threads=parameters.get("num_threads", 4),
        )
        
        # Load protein sequences
        processor = FastaProcessor()
        sequences = processor.read_file(genome_data.protein_file)
        
        # Perform BLAST search
        hits = self.search_against_reference(sequences, reference_db, blast_params)
        
        # Filter hits
        filtered_hits = self.filter_hits_by_quality(
            hits,
            min_identity=parameters.get("min_identity", 30.0),
            min_coverage=parameters.get("min_coverage", 50.0),
            max_evalue=blast_params.evalue,
        )
        
        # Cluster hits
        hit_clusters = self.cluster_hits_by_region(
            filtered_hits,
            max_distance=parameters.get("max_distance", 100000),
        )
        
        # Create analysis result
        result = AnalysisResult(
            analysis_id=f"{self.session_id}_blast",
            analysis_type=AnalysisType.BLAST_SEARCH,
            input_data={
                "genome_file": str(genome_data.genome_file),
                "protein_file": str(genome_data.protein_file),
                "organism": genome_data.organism_name,
                "reference_database": reference_db,
            },
            results={
                "blast_hits": [hit.to_dict() for hit in filtered_hits],
                "hit_clusters": [[hit.to_dict() for hit in cluster] for cluster in hit_clusters],
                "total_hits": len(hits),
                "filtered_hits": len(filtered_hits),
                "cluster_count": len(hit_clusters),
            },
            parameters=parameters or {},
            status="success",
        )
        
        self.logger.info(
            f"BLAST analysis complete: {len(filtered_hits)} hits in "
            f"{len(hit_clusters)} clusters"
        )
        
        return result
    
    def process(self, data: GenomeData, **kwargs) -> AnalysisResult:
        """Process genome data with BLAST analysis."""
        return self.analyze(data, kwargs.get('parameters'))
    
    def cleanup(self) -> None:
        """Clean up resources."""
        if self.search_engine:
            self.search_engine.cleanup()
        super().cleanup()