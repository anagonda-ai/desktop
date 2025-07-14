"""
Industry-Level Sliding Window BLAST Analysis Pipeline.

This module provides object-oriented sliding window analysis with BLAST scoring,
comprehensive error handling, parallel processing, and performance optimization.
"""

import os
import subprocess
import tempfile
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import shutil

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ...core.base_classes import (
    BaseProcessor,
    BioinformaticsConfig,
    ProcessingResult,
    SequenceProcessor,
    ParallelProcessor
)


@dataclass
class SlidingWindowConfig:
    """Configuration for sliding window analysis."""
    
    window_size: int = 10
    score_threshold: float = 0.8
    min_genes_above_threshold: int = 3
    blast_evalue: float = 1e-5
    max_target_seqs: int = 10
    blast_program: str = "blastp"  # or "blastn"
    overlap_size: int = 0  # Number of genes to overlap between windows
    
    def __post_init__(self):
        """Validate configuration."""
        if self.window_size <= 0:
            raise ValueError("Window size must be positive")
        if not 0 <= self.score_threshold <= 1:
            raise ValueError("Score threshold must be between 0 and 1")
        if self.min_genes_above_threshold <= 0:
            raise ValueError("Min genes above threshold must be positive")


@dataclass
class WindowInfo:
    """Information about a sliding window."""
    
    window_id: int
    genes: List[Dict[str, Any]] = field(default_factory=list)
    scores: List[float] = field(default_factory=list)
    sequence: str = ""
    fasta_path: Optional[Path] = None
    blast_results: Dict[str, Path] = field(default_factory=dict)
    
    @property
    def gene_count(self) -> int:
        """Get number of genes in window."""
        return len(self.genes)
    
    @property
    def high_score_count(self) -> int:
        """Get number of genes above threshold."""
        return len([s for s in self.scores if s > 0.8])  # TODO: make configurable
    
    @property
    def mean_score(self) -> float:
        """Get mean score for window."""
        return sum(self.scores) / len(self.scores) if self.scores else 0.0


@dataclass
class BlastHit:
    """BLAST hit information."""
    
    query_id: str
    subject_id: str
    identity: float
    alignment_length: int
    mismatches: int
    gap_opens: int
    q_start: int
    q_end: int
    s_start: int
    s_end: int
    evalue: float
    bit_score: float


class FastaFileManager:
    """Manager for FASTA file operations."""
    
    def __init__(self, logger):
        """Initialize FASTA file manager."""
        self.logger = logger
    
    def read_fasta_to_dataframe(self, fasta_file: Path) -> pd.DataFrame:
        """
        Read FASTA file and return DataFrame.
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            DataFrame with qseqid and sequence columns
        """
        try:
            genes = []
            
            for record in SeqIO.parse(fasta_file, "fasta"):
                genes.append({
                    "qseqid": record.id,
                    "sequence": str(record.seq)
                })
            
            df = pd.DataFrame(genes)
            
            self.logger.info(f"Loaded FASTA file", 
                           path=str(fasta_file),
                           gene_count=len(genes))
            
            return df
            
        except Exception as e:
            self.logger.error(f"Failed to read FASTA file", 
                            path=str(fasta_file),
                            error=str(e))
            raise
    
    def split_mgc_fasta_by_id(self, mgc_fasta: Path, output_dir: Path) -> List[str]:
        """
        Split MGC FASTA file by MGC IDs.
        
        Args:
            mgc_fasta: Path to MGC FASTA file
            output_dir: Output directory for split files
            
        Returns:
            List of MGC IDs
        """
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            mgc_groups = {}
            
            for record in SeqIO.parse(mgc_fasta, "fasta"):
                # Extract MGC ID (adjust based on header format)
                mgc_id = record.id.split("|")[0]
                
                if mgc_id not in mgc_groups:
                    mgc_groups[mgc_id] = []
                mgc_groups[mgc_id].append(record)
            
            # Write separate FASTA files
            for mgc_id, sequences in mgc_groups.items():
                mgc_fasta_path = output_dir / f"{mgc_id}.fasta"
                
                with open(mgc_fasta_path, "w") as out_fasta:
                    SeqIO.write(sequences, out_fasta, "fasta")
                
                self.logger.debug(f"Created MGC FASTA", 
                                mgc_id=mgc_id,
                                sequence_count=len(sequences),
                                path=str(mgc_fasta_path))
            
            self.logger.info(f"Split MGC FASTA", 
                           input_file=str(mgc_fasta),
                           output_dir=str(output_dir),
                           mgc_count=len(mgc_groups))
            
            return list(mgc_groups.keys())
            
        except Exception as e:
            self.logger.error(f"Failed to split MGC FASTA", 
                            path=str(mgc_fasta),
                            error=str(e))
            raise


class BlastRunner:
    """Professional BLAST execution and result parsing."""
    
    def __init__(self, config: SlidingWindowConfig, logger):
        """Initialize BLAST runner."""
        self.config = config
        self.logger = logger
    
    def run_blast(self, 
                  query_fasta: Path, 
                  subject_fasta: Path, 
                  output_file: Path) -> bool:
        """
        Run BLAST analysis.
        
        Args:
            query_fasta: Query FASTA file
            subject_fasta: Subject FASTA file
            output_file: Output file path
            
        Returns:
            True if successful
        """
        try:
            # Build BLAST command
            cmd = [
                self.config.blast_program,
                "-query", str(query_fasta),
                "-subject", str(subject_fasta),
                "-out", str(output_file),
                "-outfmt", "6",
                "-evalue", str(self.config.blast_evalue),
                "-max_target_seqs", str(self.config.max_target_seqs)
            ]
            
            # Run BLAST
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            if result.returncode != 0:
                self.logger.error(f"BLAST failed", 
                                command=" ".join(cmd),
                                stderr=result.stderr)
                return False
            
            self.logger.debug(f"BLAST completed successfully", 
                            query=str(query_fasta),
                            subject=str(subject_fasta),
                            output=str(output_file))
            
            return True
            
        except subprocess.TimeoutExpired:
            self.logger.error(f"BLAST timeout", 
                            query=str(query_fasta),
                            subject=str(subject_fasta))
            return False
        except Exception as e:
            self.logger.error(f"BLAST execution failed", 
                            query=str(query_fasta),
                            subject=str(subject_fasta),
                            error=str(e))
            return False
    
    def parse_blast_results(self, blast_output: Path) -> List[BlastHit]:
        """
        Parse BLAST output file.
        
        Args:
            blast_output: Path to BLAST output file
            
        Returns:
            List of BlastHit objects
        """
        try:
            hits = []
            
            if not blast_output.exists() or blast_output.stat().st_size == 0:
                return hits
            
            with open(blast_output, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 12:
                        hits.append(BlastHit(
                            query_id=parts[0],
                            subject_id=parts[1],
                            identity=float(parts[2]),
                            alignment_length=int(parts[3]),
                            mismatches=int(parts[4]),
                            gap_opens=int(parts[5]),
                            q_start=int(parts[6]),
                            q_end=int(parts[7]),
                            s_start=int(parts[8]),
                            s_end=int(parts[9]),
                            evalue=float(parts[10]),
                            bit_score=float(parts[11])
                        ))
            
            self.logger.debug(f"Parsed BLAST results", 
                            file=str(blast_output),
                            hit_count=len(hits))
            
            return hits
            
        except Exception as e:
            self.logger.error(f"Failed to parse BLAST results", 
                            file=str(blast_output),
                            error=str(e))
            return []


class SlidingWindowGenerator:
    """Generator for sliding windows with gene selection."""
    
    def __init__(self, config: SlidingWindowConfig, logger):
        """Initialize sliding window generator."""
        self.config = config
        self.logger = logger
    
    def create_windows(self, 
                      genes_df: pd.DataFrame, 
                      blast_scores_df: pd.DataFrame) -> List[WindowInfo]:
        """
        Create sliding windows with gene selection based on BLAST scores.
        
        Args:
            genes_df: DataFrame with gene information
            blast_scores_df: DataFrame with BLAST scores
            
        Returns:
            List of WindowInfo objects
        """
        try:
            # Merge genes with scores
            merged_df = pd.merge(
                genes_df, 
                blast_scores_df[['qseqid', 'normalized_composite_score']], 
                on='qseqid', 
                how='left'
            )
            merged_df['normalized_composite_score'] = merged_df['normalized_composite_score'].fillna(0)
            
            windows = []
            current_window_genes = []
            seen_genes = set()
            window_id = 0
            
            for _, row in merged_df.iterrows():
                # Extract base gene ID (remove transcript suffix)
                gene_id = row['qseqid'].split('.')[0]
                
                if gene_id not in seen_genes:
                    # Find best transcript for this gene
                    same_gene_transcripts = merged_df[
                        merged_df['qseqid'].str.startswith(gene_id)
                    ]
                    
                    if not same_gene_transcripts.empty:
                        best_idx = same_gene_transcripts['normalized_composite_score'].idxmax()
                        best_transcript = same_gene_transcripts.loc[best_idx]
                        
                        current_window_genes.append({
                            'qseqid': best_transcript['qseqid'],
                            'sequence': best_transcript['sequence']
                        })
                        
                        seen_genes.add(gene_id)
                
                # Create window when size reached
                if len(current_window_genes) == self.config.window_size:
                    window = self._create_window(window_id, current_window_genes, merged_df)
                    windows.append(window)
                    
                    # Reset for next window (with overlap if configured)
                    if self.config.overlap_size > 0:
                        overlap_genes = current_window_genes[-self.config.overlap_size:]
                        overlap_gene_ids = {g['qseqid'].split('.')[0] for g in overlap_genes}
                        current_window_genes = overlap_genes
                        seen_genes = overlap_gene_ids
                    else:
                        current_window_genes = []
                        seen_genes = set()
                    
                    window_id += 1
            
            # Add final partial window if exists
            if current_window_genes:
                window = self._create_window(window_id, current_window_genes, merged_df)
                windows.append(window)
            
            self.logger.info(f"Created sliding windows", 
                           window_count=len(windows),
                           window_size=self.config.window_size)
            
            return windows
            
        except Exception as e:
            self.logger.error(f"Failed to create sliding windows", error=str(e))
            raise
    
    def _create_window(self, 
                      window_id: int, 
                      window_genes: List[Dict[str, Any]], 
                      scores_df: pd.DataFrame) -> WindowInfo:
        """Create WindowInfo object from genes."""
        # Extract scores for window genes
        gene_ids = [g['qseqid'] for g in window_genes]
        window_scores_df = scores_df[scores_df['qseqid'].isin(gene_ids)]
        scores = window_scores_df['normalized_composite_score'].tolist()
        
        # Concatenate sequences
        sequences = [g['sequence'] for g in window_genes]
        combined_sequence = "".join(sequences)
        
        return WindowInfo(
            window_id=window_id,
            genes=window_genes,
            scores=scores,
            sequence=combined_sequence
        )


class SlidingWindowBLASTAnalyzer(SequenceProcessor):
    """
    Professional sliding window BLAST analyzer with comprehensive functionality.
    """
    
    def __init__(self, 
                 config: BioinformaticsConfig,
                 window_config: SlidingWindowConfig):
        """
        Initialize sliding window BLAST analyzer.
        
        Args:
            config: Base processing configuration
            window_config: Sliding window configuration
        """
        super().__init__(config, name="SlidingWindowBLASTAnalyzer")
        
        self.window_config = window_config
        self.fasta_manager = FastaFileManager(self.logger)
        self.blast_runner = BlastRunner(window_config, self.logger)
        self.window_generator = SlidingWindowGenerator(window_config, self.logger)
        
        # Statistics
        self.analysis_stats = {
            "windows_created": 0,
            "windows_filtered": 0,
            "blast_jobs_run": 0,
            "blast_jobs_successful": 0,
            "mgc_comparisons": 0,
        }
    
    def process(self, input_data: Dict[str, Path], **kwargs) -> ProcessingResult:
        """
        Process sliding window BLAST analysis.
        
        Args:
            input_data: Dictionary with required file paths:
                - fasta_file: Path to input FASTA file
                - blast_file: Path to BLAST scores file
                - mgc_fasta: Path to MGC FASTA file
            **kwargs: Additional processing parameters
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            self.logger.info("Starting sliding window BLAST analysis")
            
            # Validate inputs
            required_files = ['fasta_file', 'blast_file', 'mgc_fasta']
            for file_key in required_files:
                if file_key not in input_data:
                    result.add_error(f"Missing required input: {file_key}")
                    return result
                
                if not input_data[file_key].exists():
                    result.add_error(f"Input file does not exist: {input_data[file_key]}")
                    return result
            
            # Create output directory
            output_dir = self.config.output_dir / "sliding_window_analysis"
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Step 1: Load data
            genes_df = self.fasta_manager.read_fasta_to_dataframe(input_data['fasta_file'])
            blast_scores_df = pd.read_csv(input_data['blast_file'], sep='\t')
            
            # Step 2: Create sliding windows
            windows = self.window_generator.create_windows(genes_df, blast_scores_df)
            self.analysis_stats["windows_created"] = len(windows)
            
            # Step 3: Filter windows
            filtered_windows = self._filter_windows(windows)
            self.analysis_stats["windows_filtered"] = len(filtered_windows)
            
            # Step 4: Split MGC FASTA
            mgc_dir = output_dir / "MGC_FASTAs"
            mgc_ids = self.fasta_manager.split_mgc_fasta_by_id(input_data['mgc_fasta'], mgc_dir)
            
            # Step 5: Run BLAST analysis
            blast_results = self._run_blast_analysis(filtered_windows, mgc_dir, output_dir)
            
            # Step 6: Generate results summary
            summary = self._generate_results_summary(filtered_windows, blast_results, output_dir)
            
            result.metadata.update({
                "windows_created": len(windows),
                "windows_filtered": len(filtered_windows),
                "mgc_count": len(mgc_ids),
                "blast_comparisons": self.analysis_stats["mgc_comparisons"],
                "output_directory": str(output_dir),
                "results_summary": summary,
                "analysis_stats": self.analysis_stats
            })
            
            self.logger.info(f"Sliding window BLAST analysis completed",
                           windows_created=len(windows),
                           windows_filtered=len(filtered_windows),
                           mgc_count=len(mgc_ids))
            
        except Exception as e:
            result.add_error(f"Analysis failed: {str(e)}")
        
        return result
    
    def _filter_windows(self, windows: List[WindowInfo]) -> List[WindowInfo]:
        """
        Filter windows based on score criteria.
        
        Args:
            windows: List of windows to filter
            
        Returns:
            List of filtered windows
        """
        filtered = []
        
        for window in windows:
            genes_above_threshold = sum(
                1 for score in window.scores 
                if score > self.window_config.score_threshold
            )
            
            if genes_above_threshold >= self.window_config.min_genes_above_threshold:
                filtered.append(window)
        
        self.logger.info(f"Filtered windows", 
                       original=len(windows),
                       filtered=len(filtered),
                       threshold=self.window_config.score_threshold)
        
        return filtered
    
    def _run_blast_analysis(self, 
                           windows: List[WindowInfo], 
                           mgc_dir: Path, 
                           output_dir: Path) -> Dict[int, Dict[str, Path]]:
        """
        Run BLAST analysis for all windows against all MGCs.
        
        Args:
            windows: List of filtered windows
            mgc_dir: Directory containing MGC FASTA files
            output_dir: Output directory
            
        Returns:
            Dictionary mapping window IDs to BLAST result files
        """
        blast_results_dir = output_dir / "BLAST_results"
        blast_results_dir.mkdir(parents=True, exist_ok=True)
        
        mgc_files = list(mgc_dir.glob("*.fasta"))
        all_results = {}
        
        for window in windows:
            window_dir = blast_results_dir / f"window_{window.window_id}"
            window_dir.mkdir(parents=True, exist_ok=True)
            
            # Create window FASTA file
            window_fasta = window_dir / f"window_{window.window_id}_query.fasta"
            self._create_window_fasta(window, window_fasta)
            
            window_results = {}
            
            # BLAST against each MGC
            for mgc_fasta in mgc_files:
                mgc_id = mgc_fasta.stem
                blast_output = window_dir / f"window_{window.window_id}_vs_{mgc_id}_blast.txt"
                
                if self.blast_runner.run_blast(window_fasta, mgc_fasta, blast_output):
                    window_results[mgc_id] = blast_output
                    self.analysis_stats["blast_jobs_successful"] += 1
                
                self.analysis_stats["blast_jobs_run"] += 1
                self.analysis_stats["mgc_comparisons"] += 1
            
            all_results[window.window_id] = window_results
            
            self.logger.debug(f"Completed BLAST for window", 
                            window_id=window.window_id,
                            mgc_count=len(window_results))
        
        return all_results
    
    def _create_window_fasta(self, window: WindowInfo, output_path: Path) -> None:
        """Create FASTA file for window."""
        try:
            with open(output_path, "w") as f:
                f.write(f">window_{window.window_id}_query\n")
                f.write(f"{window.sequence}\n")
            
            window.fasta_path = output_path
            
        except Exception as e:
            self.logger.error(f"Failed to create window FASTA", 
                            window_id=window.window_id,
                            error=str(e))
            raise
    
    def _generate_results_summary(self, 
                                 windows: List[WindowInfo], 
                                 blast_results: Dict[int, Dict[str, Path]], 
                                 output_dir: Path) -> Path:
        """Generate comprehensive results summary."""
        try:
            summary_data = []
            
            for window in windows:
                window_id = window.window_id
                window_blast_results = blast_results.get(window_id, {})
                
                for mgc_id, result_file in window_blast_results.items():
                    # Parse BLAST results
                    hits = self.blast_runner.parse_blast_results(result_file)
                    
                    summary_data.append({
                        "window_id": window_id,
                        "mgc_id": mgc_id,
                        "gene_count": window.gene_count,
                        "mean_score": window.mean_score,
                        "high_score_count": window.high_score_count,
                        "blast_hits": len(hits),
                        "best_bit_score": max([h.bit_score for h in hits], default=0),
                        "best_evalue": min([h.evalue for h in hits], default=float('inf')),
                        "result_file": str(result_file)
                    })
            
            # Save summary
            summary_file = output_dir / "analysis_summary.csv"
            summary_df = pd.DataFrame(summary_data)
            summary_df.to_csv(summary_file, index=False)
            
            self.logger.info(f"Generated results summary", 
                           file=str(summary_file),
                           records=len(summary_data))
            
            return summary_file
            
        except Exception as e:
            self.logger.error(f"Failed to generate results summary", error=str(e))
            raise
    
    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data.
        
        Args:
            input_data: Input data to validate
            
        Returns:
            True if valid
        """
        if not isinstance(input_data, dict):
            return False
        
        required_keys = ['fasta_file', 'blast_file', 'mgc_fasta']
        return all(key in input_data for key in required_keys)


def main():
    """Main entry point for sliding window BLAST analysis."""
    # File paths
    input_files = {
        'fasta_file': Path("/groups/itay_mayrose/alongonda/desktop/arabidopsis/sorted_Arabidopsis_thaliana.TAIR10.pep.all.fa"),
        'blast_file': Path("/groups/itay_mayrose/alongonda/desktop/Arabidopsis_thaliana.TAIR10.pep.ordered.fa/best_normalized_blast_scores.csv"),
        'mgc_fasta': Path("/groups/itay_mayrose/alongonda/desktop/MGCs/all_genes_from_mibig/mibig_prot_seqs_4.0.fasta")
    }
    
    # Configuration
    config = BioinformaticsConfig(
        input_dir=Path("/groups/itay_mayrose/alongonda/desktop"),
        output_dir=Path("/groups/itay_mayrose/alongonda/desktop/sliding_window_output"),
        max_workers=8,
        log_level="INFO"
    )
    
    window_config = SlidingWindowConfig(
        window_size=10,
        score_threshold=0.8,
        min_genes_above_threshold=3,
        blast_evalue=1e-5,
        blast_program="blastp"
    )
    
    # Initialize analyzer
    try:
        analyzer = SlidingWindowBLASTAnalyzer(
            config=config,
            window_config=window_config
        )
        
        # Run analysis
        result = analyzer.run(input_files)
        
        # Display results
        if result.success:
            print("✅ Sliding window BLAST analysis completed successfully!")
            print(f"   Windows created: {result.metadata.get('windows_created', 0)}")
            print(f"   Windows filtered: {result.metadata.get('windows_filtered', 0)}")
            print(f"   MGC comparisons: {result.metadata.get('blast_comparisons', 0)}")
            print(f"   Processing time: {result.processing_time:.2f} seconds")
            print(f"   Output directory: {result.metadata.get('output_directory')}")
        else:
            print("❌ Sliding window BLAST analysis failed!")
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