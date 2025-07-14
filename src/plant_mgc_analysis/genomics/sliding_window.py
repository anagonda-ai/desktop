"""
Sliding Window Analysis for Plant MGC Analysis Pipeline.

This module provides object-oriented sliding window analysis for identifying
metabolic gene clusters in plant genomes with proper statistical validation.
"""

import statistics
from typing import Dict, List, Optional, Tuple, Iterator, Set
from dataclasses import dataclass, field
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np
from scipy import stats
from loguru import logger

from ..core.base import AnalysisEngine, StatisticalAnalyzer
from ..core.types import (
    GenomeData,
    GeneInfo, 
    MGCCandidate,
    AnalysisResult,
    AnalysisType,
    MGCStatus,
)
from ..core.exceptions import AnalysisError, ValidationError
from ..utils.file_operations import CsvProcessor, GffProcessor
from ..utils.validation import validate_analysis_parameters


@dataclass
class SlidingWindow:
    """Represents a genomic sliding window."""
    
    chromosome: str
    start: int
    end: int
    genes: List[GeneInfo] = field(default_factory=list)
    metabolic_genes: List[GeneInfo] = field(default_factory=list)
    
    @property
    def size(self) -> int:
        """Return window size in base pairs."""
        return self.end - self.start
    
    @property
    def gene_count(self) -> int:
        """Return total number of genes in window."""
        return len(self.genes)
    
    @property
    def metabolic_gene_count(self) -> int:
        """Return number of metabolic genes in window."""
        return len(self.metabolic_genes)
    
    @property
    def gene_density(self) -> float:
        """Return gene density (genes per kb)."""
        if self.size == 0:
            return 0.0
        return (self.gene_count / self.size) * 1000
    
    @property
    def metabolic_enrichment(self) -> float:
        """Return metabolic gene enrichment ratio."""
        if self.gene_count == 0:
            return 0.0
        return self.metabolic_gene_count / self.gene_count
    
    def get_ec_families(self) -> Set[str]:
        """Get unique EC number families in this window."""
        families = set()
        for gene in self.metabolic_genes:
            if hasattr(gene, 'ec_number') and gene.ec_number:
                family = gene.ec_number.split('.')[0]
                families.add(family)
        return families
    
    def get_ec_subfamilies(self) -> Set[str]:
        """Get unique EC number subfamilies in this window."""
        subfamilies = set()
        for gene in self.metabolic_genes:
            if hasattr(gene, 'ec_number') and gene.ec_number:
                subfamily = '.'.join(gene.ec_number.split('.')[:2])
                subfamilies.add(subfamily)
        return subfamilies


@dataclass
class EnrichmentStatistics:
    """Statistical measures for metabolic gene enrichment."""
    
    enrichment_ratio: float
    p_value: float
    z_score: float
    confidence_interval: Tuple[float, float]
    background_rate: float
    observed_count: int
    expected_count: float
    window_gene_count: int
    
    @property
    def is_significant(self) -> bool:
        """Check if enrichment is statistically significant."""
        return self.p_value < 0.05
    
    @property
    def is_enriched(self) -> bool:
        """Check if window is enriched for metabolic genes."""
        return self.enrichment_ratio > 1.0 and self.is_significant


class MetabolicGeneDatabase:
    """Database interface for metabolic gene annotations."""
    
    def __init__(self):
        """Initialize metabolic gene database."""
        self._gene_annotations: Dict[str, Dict[str, str]] = {}
        self._ec_to_pathways: Dict[str, List[str]] = {}
        self._pathway_names: Dict[str, str] = {}
    
    def load_annotations(self, annotation_file: Path) -> None:
        """
        Load gene annotations from file.
        
        Args:
            annotation_file: Path to annotation file (CSV format)
        """
        processor = CsvProcessor()
        df = processor.read_file(annotation_file)
        
        required_columns = ['gene_id', 'ec_number']
        if not all(col in df.columns for col in required_columns):
            raise ValidationError(
                f"Annotation file must contain columns: {required_columns}",
                field_name="file_columns",
                field_value=list(df.columns)
            )
        
        for _, row in df.iterrows():
            gene_id = str(row['gene_id']).lower()
            self._gene_annotations[gene_id] = {
                'ec_number': row.get('ec_number', ''),
                'pathway': row.get('pathway', ''),
                'description': row.get('description', ''),
                'metabolic_role': row.get('metabolic_role', ''),
            }
        
        logger.info(f"Loaded annotations for {len(self._gene_annotations)} genes")
    
    def is_metabolic_gene(self, gene_id: str) -> bool:
        """Check if gene is metabolic."""
        return gene_id.lower() in self._gene_annotations
    
    def get_annotation(self, gene_id: str) -> Optional[Dict[str, str]]:
        """Get annotation for gene."""
        return self._gene_annotations.get(gene_id.lower())
    
    def get_ec_number(self, gene_id: str) -> Optional[str]:
        """Get EC number for gene."""
        annotation = self.get_annotation(gene_id)
        return annotation.get('ec_number') if annotation else None
    
    def get_pathway(self, gene_id: str) -> Optional[str]:
        """Get pathway for gene."""
        annotation = self.get_annotation(gene_id)
        return annotation.get('pathway') if annotation else None


class SlidingWindowAnalyzer(AnalysisEngine, StatisticalAnalyzer):
    """
    Sliding window analyzer for metabolic gene cluster detection.
    
    This class implements sophisticated sliding window analysis with statistical
    validation to identify regions enriched for metabolic genes.
    """
    
    def __init__(
        self,
        window_size: int = 50000,
        step_size: int = 10000,
        min_genes_per_window: int = 3,
        min_metabolic_genes: int = 2,
        min_ec_families: int = 2,
        min_ec_subfamilies: int = 2,
        **kwargs
    ):
        """
        Initialize sliding window analyzer.
        
        Args:
            window_size: Size of sliding window in base pairs
            step_size: Step size for sliding window
            min_genes_per_window: Minimum genes required per window
            min_metabolic_genes: Minimum metabolic genes for candidacy
            min_ec_families: Minimum EC number families required
            min_ec_subfamilies: Minimum EC number subfamilies required
            **kwargs: Additional arguments for parent classes
        """
        super().__init__(
            analysis_type=AnalysisType.SLIDING_WINDOW,
            **kwargs
        )
        
        self.window_size = window_size
        self.step_size = step_size
        self.min_genes_per_window = min_genes_per_window
        self.min_metabolic_genes = min_metabolic_genes
        self.min_ec_families = min_ec_families
        self.min_ec_subfamilies = min_ec_subfamilies
        
        self.metabolic_db = MetabolicGeneDatabase()
        self._genome_genes: Dict[str, List[GeneInfo]] = {}
        self._background_metabolic_rate: float = 0.0
    
    def validate_input(self, genome_data: GenomeData) -> None:
        """Validate input genome data."""
        if not genome_data.annotation_file:
            raise ValidationError(
                "Annotation file required for sliding window analysis",
                field_name="annotation_file",
                field_value=None
            )
        
        if not genome_data.annotation_file.exists():
            raise ValidationError(
                f"Annotation file not found: {genome_data.annotation_file}",
                field_name="annotation_file",
                field_value=str(genome_data.annotation_file)
            )
    
    def validate_parameters(self, parameters: Dict[str, str]) -> None:
        """Validate analysis parameters."""
        super().validate_parameters(parameters)
        
        validated = validate_analysis_parameters("sliding_window", parameters)
        
        # Update instance parameters
        self.window_size = validated.get("window_size", self.window_size)
        self.step_size = validated.get("step_size", self.step_size)
        self.min_genes_per_window = validated.get("min_genes", self.min_genes_per_window)
    
    def load_metabolic_annotations(self, annotation_file: Path) -> None:
        """Load metabolic gene annotations."""
        self.metabolic_db.load_annotations(annotation_file)
    
    def load_genome_annotation(self, genome_data: GenomeData) -> None:
        """Load genome annotation data."""
        processor = GffProcessor()
        genes = processor.read_file(genome_data.annotation_file)
        
        # Group genes by chromosome
        for gene in genes:
            if gene.chromosome not in self._genome_genes:
                self._genome_genes[gene.chromosome] = []
            self._genome_genes[gene.chromosome].append(gene)
        
        # Sort genes by position on each chromosome
        for chromosome in self._genome_genes:
            self._genome_genes[chromosome].sort(key=lambda g: g.start)
        
        # Calculate background metabolic gene rate
        total_genes = sum(len(genes) for genes in self._genome_genes.values())
        metabolic_genes = sum(
            1 for genes in self._genome_genes.values() 
            for gene in genes 
            if self.metabolic_db.is_metabolic_gene(gene.gene_id)
        )
        
        self._background_metabolic_rate = metabolic_genes / total_genes if total_genes > 0 else 0.0
        
        logger.info(
            f"Loaded {total_genes} genes across {len(self._genome_genes)} chromosomes. "
            f"Background metabolic rate: {self._background_metabolic_rate:.3f}"
        )
    
    def generate_windows(self, chromosome: str) -> Iterator[SlidingWindow]:
        """
        Generate sliding windows for a chromosome.
        
        Args:
            chromosome: Chromosome name
            
        Yields:
            SlidingWindow objects
        """
        if chromosome not in self._genome_genes:
            return
        
        genes = self._genome_genes[chromosome]
        if not genes:
            return
        
        # Get chromosome boundaries
        chr_start = min(gene.start for gene in genes)
        chr_end = max(gene.end for gene in genes)
        
        # Generate windows
        window_start = chr_start
        while window_start < chr_end:
            window_end = window_start + self.window_size
            
            # Find genes in this window
            window_genes = [
                gene for gene in genes
                if not (gene.end < window_start or gene.start > window_end)
            ]
            
            # Identify metabolic genes
            metabolic_genes = []
            for gene in window_genes:
                if self.metabolic_db.is_metabolic_gene(gene.gene_id):
                    # Add EC number to gene object
                    ec_number = self.metabolic_db.get_ec_number(gene.gene_id)
                    if ec_number:
                        gene.ec_number = ec_number
                        gene.metabolic_role = 'metabolic'
                        metabolic_genes.append(gene)
            
            window = SlidingWindow(
                chromosome=chromosome,
                start=window_start,
                end=window_end,
                genes=window_genes,
                metabolic_genes=metabolic_genes,
            )
            
            yield window
            window_start += self.step_size
    
    def calculate_enrichment_statistics(self, window: SlidingWindow) -> EnrichmentStatistics:
        """
        Calculate enrichment statistics for a window.
        
        Args:
            window: Sliding window to analyze
            
        Returns:
            Enrichment statistics
        """
        observed_count = window.metabolic_gene_count
        total_count = window.gene_count
        expected_count = total_count * self._background_metabolic_rate
        
        # Calculate enrichment ratio
        if expected_count > 0:
            enrichment_ratio = observed_count / expected_count
        else:
            enrichment_ratio = float('inf') if observed_count > 0 else 1.0
        
        # Perform binomial test
        if total_count > 0 and self._background_metabolic_rate > 0:
            p_value = stats.binom_test(
                observed_count,
                total_count,
                self._background_metabolic_rate,
                alternative='greater'
            )
        else:
            p_value = 1.0
        
        # Calculate z-score
        if total_count > 0 and self._background_metabolic_rate > 0:
            variance = total_count * self._background_metabolic_rate * (1 - self._background_metabolic_rate)
            if variance > 0:
                z_score = (observed_count - expected_count) / np.sqrt(variance)
            else:
                z_score = 0.0
        else:
            z_score = 0.0
        
        # Calculate confidence interval (95%)
        if total_count > 0:
            conf_int = stats.binom.interval(
                0.95, total_count, self._background_metabolic_rate
            )
        else:
            conf_int = (0.0, 0.0)
        
        return EnrichmentStatistics(
            enrichment_ratio=enrichment_ratio,
            p_value=p_value,
            z_score=z_score,
            confidence_interval=conf_int,
            background_rate=self._background_metabolic_rate,
            observed_count=observed_count,
            expected_count=expected_count,
            window_gene_count=total_count,
        )
    
    def filter_candidate_windows(self, window: SlidingWindow, stats: EnrichmentStatistics) -> bool:
        """
        Filter windows for MGC candidacy.
        
        Args:
            window: Sliding window
            stats: Enrichment statistics
            
        Returns:
            True if window passes filtering criteria
        """
        # Basic filters
        if window.gene_count < self.min_genes_per_window:
            return False
        
        if window.metabolic_gene_count < self.min_metabolic_genes:
            return False
        
        if not stats.is_significant:
            return False
        
        # EC number diversity filters
        ec_families = window.get_ec_families()
        ec_subfamilies = window.get_ec_subfamilies()
        
        if len(ec_families) < self.min_ec_families:
            return False
        
        # Apply subfamily filtering logic from original code
        if len(ec_subfamilies) >= 3:
            return True
        elif len(ec_subfamilies) == 2:
            # Check if subfamilies span multiple families
            subfamily_families = defaultdict(set)
            for gene in window.metabolic_genes:
                if hasattr(gene, 'ec_number') and gene.ec_number:
                    family = gene.ec_number.split('.')[0]
                    subfamily = '.'.join(gene.ec_number.split('.')[:2])
                    subfamily_families[subfamily].add(family)
            
            return all(len(families) > 1 for families in subfamily_families.values())
        
        return False
    
    def create_mgc_candidate(
        self, 
        window: SlidingWindow, 
        stats: EnrichmentStatistics,
        candidate_id: str,
        organism: str
    ) -> MGCCandidate:
        """
        Create MGC candidate from window.
        
        Args:
            window: Sliding window
            stats: Enrichment statistics
            candidate_id: Unique candidate identifier
            organism: Organism name
            
        Returns:
            MGC candidate object
        """
        # Calculate actual cluster boundaries from genes
        if window.metabolic_genes:
            actual_start = min(gene.start for gene in window.metabolic_genes)
            actual_end = max(gene.end for gene in window.metabolic_genes)
        else:
            actual_start = window.start
            actual_end = window.end
        
        # Determine pathway information
        pathways = set()
        for gene in window.metabolic_genes:
            pathway = self.metabolic_db.get_pathway(gene.gene_id)
            if pathway:
                pathways.add(pathway)
        
        primary_pathway = list(pathways)[0] if pathways else None
        
        return MGCCandidate(
            cluster_id=candidate_id,
            chromosome=window.chromosome,
            start=actual_start,
            end=actual_end,
            genes=window.genes,
            organism=organism,
            analysis_method="sliding_window",
            confidence_score=1.0 - stats.p_value,
            status=MGCStatus.CANDIDATE,
            pathway_id=primary_pathway,
            pathway_name=primary_pathway,
            metabolic_class="mixed" if len(pathways) > 1 else "single_pathway",
            validation_data={
                "enrichment_ratio": stats.enrichment_ratio,
                "p_value": stats.p_value,
                "z_score": stats.z_score,
                "metabolic_gene_count": window.metabolic_gene_count,
                "total_gene_count": window.gene_count,
                "ec_families": len(window.get_ec_families()),
                "ec_subfamilies": len(window.get_ec_subfamilies()),
                "window_size": window.size,
                "gene_density": window.gene_density,
            }
        )
    
    def analyze(
        self,
        genome_data: GenomeData,
        parameters: Optional[Dict[str, str]] = None
    ) -> AnalysisResult:
        """
        Perform sliding window analysis.
        
        Args:
            genome_data: Genome data to analyze
            parameters: Optional analysis parameters
            
        Returns:
            Analysis result
        """
        if parameters:
            self.validate_parameters(parameters)
        
        # Load genome annotation
        self.load_genome_annotation(genome_data)
        
        # Load metabolic annotations if provided
        if 'metabolic_annotation_file' in (parameters or {}):
            annotation_file = Path(parameters['metabolic_annotation_file'])
            self.load_metabolic_annotations(annotation_file)
        
        candidates = []
        all_windows = []
        candidate_count = 0
        
        # Process each chromosome
        for chromosome in self._genome_genes:
            self.logger.info(f"Processing chromosome {chromosome}")
            
            for window in self.generate_windows(chromosome):
                # Calculate enrichment statistics
                stats = self.calculate_enrichment_statistics(window)
                
                # Store window data for analysis
                all_windows.append({
                    "chromosome": window.chromosome,
                    "start": window.start,
                    "end": window.end,
                    "gene_count": window.gene_count,
                    "metabolic_gene_count": window.metabolic_gene_count,
                    "enrichment_ratio": stats.enrichment_ratio,
                    "p_value": stats.p_value,
                    "z_score": stats.z_score,
                })
                
                # Check if window passes filtering
                if self.filter_candidate_windows(window, stats):
                    candidate_count += 1
                    candidate_id = f"{genome_data.organism_name}_SW_{candidate_count:04d}"
                    
                    candidate = self.create_mgc_candidate(
                        window, stats, candidate_id, genome_data.organism_name
                    )
                    candidates.append(candidate)
        
        # Apply multiple testing correction
        if all_windows:
            p_values = [w["p_value"] for w in all_windows]
            corrected_p_values = self.apply_multiple_testing_correction(p_values)
            
            for window, corrected_p in zip(all_windows, corrected_p_values):
                window["corrected_p_value"] = corrected_p
        
        # Create analysis result
        result = AnalysisResult(
            analysis_id=f"{self.session_id}_sliding_window",
            analysis_type=AnalysisType.SLIDING_WINDOW,
            input_data={
                "genome_file": str(genome_data.genome_file),
                "annotation_file": str(genome_data.annotation_file),
                "organism": genome_data.organism_name,
                "window_size": self.window_size,
                "step_size": self.step_size,
            },
            results={
                "mgc_candidates": candidates,
                "all_windows": all_windows,
                "total_windows": len(all_windows),
                "candidate_count": len(candidates),
                "background_metabolic_rate": self._background_metabolic_rate,
            },
            parameters=parameters or {},
            status="success",
        )
        
        self.logger.info(
            f"Sliding window analysis complete: {len(candidates)} candidates "
            f"from {len(all_windows)} windows"
        )
        
        return result
    
    def process(self, data: GenomeData, **kwargs) -> AnalysisResult:
        """Process genome data with sliding window analysis."""
        return self.analyze(data, kwargs.get('parameters'))
    
    def calculate_statistics(self, data: List[SlidingWindow]) -> Dict[str, float]:
        """Calculate summary statistics for windows."""
        if not data:
            return {}
        
        enrichments = [w.metabolic_enrichment for w in data]
        gene_counts = [w.gene_count for w in data]
        metabolic_counts = [w.metabolic_gene_count for w in data]
        
        return {
            "mean_enrichment": statistics.mean(enrichments),
            "median_enrichment": statistics.median(enrichments),
            "std_enrichment": statistics.stdev(enrichments) if len(enrichments) > 1 else 0.0,
            "mean_gene_count": statistics.mean(gene_counts),
            "mean_metabolic_count": statistics.mean(metabolic_counts),
            "total_windows": len(data),
        }
    
    def perform_hypothesis_test(
        self,
        sample1: List[float],
        sample2: List[float],
        test_type: str = "two_sided"
    ) -> Dict[str, float]:
        """Perform statistical test between two samples."""
        statistic, p_value = stats.mannwhitneyu(
            sample1, sample2, alternative=test_type
        )
        
        return {
            "statistic": statistic,
            "p_value": p_value,
            "test_type": test_type,
            "sample1_size": len(sample1),
            "sample2_size": len(sample2),
        }