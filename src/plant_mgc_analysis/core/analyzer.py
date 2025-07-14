"""
Main analyzer class for Plant MGC Analysis Pipeline.

This module provides the primary interface for conducting metabolic gene cluster
analysis on plant genomes.
"""

import uuid
from pathlib import Path
from typing import List, Dict, Any, Optional, Union
from datetime import datetime
from loguru import logger

from .types import (
    GenomeData,
    MGCCandidate,
    AnalysisResult,
    AnalysisType,
    PipelineConfig,
    PathLike,
)
from .exceptions import (
    AnalysisError,
    ValidationError,
    ConfigurationError,
)
from ..config.settings import get_settings
from ..utils.logging import setup_logging
from ..utils.validation import validate_input_files


class MGCAnalyzer:
    """
    Main analyzer class for plant metabolic gene cluster analysis.
    
    This class orchestrates the entire analysis pipeline, including:
    - Genome data processing
    - MGC candidate identification
    - Statistical analysis
    - Machine learning predictions
    - Results integration and reporting
    """
    
    def __init__(
        self,
        config: Optional[Dict[str, Any]] = None,
        log_level: str = "INFO",
    ):
        """
        Initialize MGC Analyzer.
        
        Args:
            config: Optional configuration dictionary
            log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        """
        self.settings = get_settings()
        if config:
            self.settings.update_config(**config)
        
        # Setup logging
        setup_logging(
            level=log_level,
            log_file=self.settings.logging.log_file,
            format_string=self.settings.logging.format,
        )
        
        # Initialize analysis state
        self.analysis_id = str(uuid.uuid4())
        self.genome_data: Optional[GenomeData] = None
        self.analysis_results: List[AnalysisResult] = []
        self.pipeline_config: Optional[PipelineConfig] = None
        
        logger.info(f"MGC Analyzer initialized with ID: {self.analysis_id}")
    
    def load_genome(
        self,
        genome_file: PathLike,
        annotation_file: Optional[PathLike] = None,
        protein_file: Optional[PathLike] = None,
        organism_name: Optional[str] = None,
        **metadata: Any,
    ) -> GenomeData:
        """
        Load genome data for analysis.
        
        Args:
            genome_file: Path to genome FASTA file
            annotation_file: Optional path to annotation file (GFF/GTF)
            protein_file: Optional path to protein FASTA file
            organism_name: Name of organism
            **metadata: Additional metadata
            
        Returns:
            GenomeData object
            
        Raises:
            ValidationError: If input files are invalid
            FileNotFoundError: If required files don't exist
        """
        logger.info(f"Loading genome data from {genome_file}")
        
        # Validate input files
        genome_path = Path(genome_file)
        annotation_path = Path(annotation_file) if annotation_file else None
        protein_path = Path(protein_file) if protein_file else None
        
        validate_input_files([genome_path])
        
        if annotation_path:
            validate_input_files([annotation_path])
        
        if protein_path:
            validate_input_files([protein_path])
        
        # Extract organism name from filename if not provided
        if not organism_name:
            organism_name = genome_path.stem
        
        # Create genome data object
        self.genome_data = GenomeData(
            organism_name=organism_name,
            genome_file=genome_path,
            annotation_file=annotation_path,
            protein_file=protein_path,
            metadata=metadata,
        )
        
        logger.info(f"Genome data loaded for {organism_name}")
        return self.genome_data
    
    def configure_pipeline(
        self,
        analysis_types: List[Union[str, AnalysisType]],
        output_dir: PathLike,
        parameters: Optional[Dict[str, Any]] = None,
        parallel_jobs: int = 1,
        use_cache: bool = True,
    ) -> PipelineConfig:
        """
        Configure analysis pipeline.
        
        Args:
            analysis_types: List of analysis types to perform
            output_dir: Output directory for results
            parameters: Optional analysis parameters
            parallel_jobs: Number of parallel jobs
            use_cache: Whether to use caching
            
        Returns:
            PipelineConfig object
            
        Raises:
            ValidationError: If configuration is invalid
        """
        logger.info("Configuring analysis pipeline")
        
        # Convert string analysis types to enum
        converted_types = []
        for analysis_type in analysis_types:
            if isinstance(analysis_type, str):
                try:
                    converted_types.append(AnalysisType(analysis_type))
                except ValueError:
                    raise ValidationError(
                        f"Invalid analysis type: {analysis_type}",
                        field_name="analysis_types",
                        field_value=analysis_type,
                    )
            else:
                converted_types.append(analysis_type)
        
        # Create pipeline configuration
        self.pipeline_config = PipelineConfig(
            pipeline_name=f"mgc_analysis_{self.analysis_id}",
            analysis_types=converted_types,
            input_files=[self.genome_data.genome_file] if self.genome_data else [],
            output_dir=Path(output_dir),
            parameters=parameters or {},
            parallel_jobs=parallel_jobs,
            use_cache=use_cache,
            cache_dir=self.settings.cache_dir,
        )
        
        # Validate configuration
        self.pipeline_config.validate()
        
        logger.info(f"Pipeline configured with {len(converted_types)} analysis types")
        return self.pipeline_config
    
    def run_analysis(
        self,
        analysis_type: Union[str, AnalysisType],
        parameters: Optional[Dict[str, Any]] = None,
    ) -> AnalysisResult:
        """
        Run a specific analysis.
        
        Args:
            analysis_type: Type of analysis to run
            parameters: Optional analysis parameters
            
        Returns:
            AnalysisResult object
            
        Raises:
            AnalysisError: If analysis fails
            ValidationError: If inputs are invalid
        """
        if self.genome_data is None:
            raise ValidationError("No genome data loaded. Call load_genome() first.")
        
        # Convert string to enum if needed
        if isinstance(analysis_type, str):
            try:
                analysis_type = AnalysisType(analysis_type)
            except ValueError:
                raise ValidationError(
                    f"Invalid analysis type: {analysis_type}",
                    field_name="analysis_type",
                    field_value=analysis_type,
                )
        
        logger.info(f"Running {analysis_type.value} analysis")
        
        # Create analysis result
        result = AnalysisResult(
            analysis_id=f"{self.analysis_id}_{analysis_type.value}",
            analysis_type=analysis_type,
            input_data={
                "genome_file": str(self.genome_data.genome_file),
                "organism": self.genome_data.organism_name,
            },
            results={},
            parameters=parameters or {},
            status="running",
        )
        
        start_time = datetime.now()
        
        try:
            # Route to appropriate analysis method
            if analysis_type == AnalysisType.SLIDING_WINDOW:
                result.results = self._run_sliding_window_analysis(parameters)
            elif analysis_type == AnalysisType.BLAST_SEARCH:
                result.results = self._run_blast_analysis(parameters)
            elif analysis_type == AnalysisType.PHYLOGENETIC:
                result.results = self._run_phylogenetic_analysis(parameters)
            elif analysis_type == AnalysisType.MACHINE_LEARNING:
                result.results = self._run_ml_analysis(parameters)
            elif analysis_type == AnalysisType.PATHWAY_ANALYSIS:
                result.results = self._run_pathway_analysis(parameters)
            elif analysis_type == AnalysisType.COMPARATIVE_GENOMICS:
                result.results = self._run_comparative_analysis(parameters)
            else:
                raise AnalysisError(f"Unsupported analysis type: {analysis_type}")
            
            # Update result status
            result.status = "success"
            result.execution_time = (datetime.now() - start_time).total_seconds()
            
            logger.info(f"Analysis {analysis_type.value} completed successfully")
            
        except Exception as e:
            result.status = "failed"
            result.error_message = str(e)
            result.execution_time = (datetime.now() - start_time).total_seconds()
            
            logger.error(f"Analysis {analysis_type.value} failed: {e}")
            raise AnalysisError(
                f"Analysis {analysis_type.value} failed: {e}",
                analysis_type=analysis_type.value,
                input_data=str(self.genome_data.genome_file),
            )
        
        # Store result
        self.analysis_results.append(result)
        return result
    
    def run_pipeline(self) -> List[AnalysisResult]:
        """
        Run the complete analysis pipeline.
        
        Returns:
            List of AnalysisResult objects
            
        Raises:
            ConfigurationError: If pipeline is not configured
            AnalysisError: If any analysis fails
        """
        if self.pipeline_config is None:
            raise ConfigurationError("Pipeline not configured. Call configure_pipeline() first.")
        
        logger.info("Starting analysis pipeline")
        
        results = []
        for analysis_type in self.pipeline_config.analysis_types:
            try:
                result = self.run_analysis(
                    analysis_type,
                    self.pipeline_config.parameters.get(analysis_type.value, {})
                )
                results.append(result)
            except Exception as e:
                logger.error(f"Pipeline failed at {analysis_type.value}: {e}")
                raise
        
        logger.info("Analysis pipeline completed")
        return results
    
    def get_mgc_candidates(self) -> List[MGCCandidate]:
        """
        Get all MGC candidates from analysis results.
        
        Returns:
            List of MGCCandidate objects
        """
        candidates = []
        for result in self.analysis_results:
            candidates.extend(result.get_mgc_candidates())
        return candidates
    
    def save_results(self, output_file: PathLike) -> None:
        """
        Save analysis results to file.
        
        Args:
            output_file: Path to output file
        """
        # Implementation will be added in utils module
        logger.info(f"Saving results to {output_file}")
        # TODO: Implement result saving
    
    def generate_report(self) -> Dict[str, Any]:
        """
        Generate analysis report.
        
        Returns:
            Dictionary containing report data
        """
        logger.info("Generating analysis report")
        
        candidates = self.get_mgc_candidates()
        
        report = {
            "analysis_id": self.analysis_id,
            "organism": self.genome_data.organism_name if self.genome_data else "Unknown",
            "analysis_types": [r.analysis_type.value for r in self.analysis_results],
            "total_candidates": len(candidates),
            "successful_analyses": len([r for r in self.analysis_results if r.success]),
            "failed_analyses": len([r for r in self.analysis_results if not r.success]),
            "execution_time": sum(r.execution_time or 0 for r in self.analysis_results),
            "created_at": datetime.now().isoformat(),
            "results_summary": {
                result.analysis_type.value: {
                    "status": result.status,
                    "candidates": len(result.get_mgc_candidates()),
                    "execution_time": result.execution_time,
                }
                for result in self.analysis_results
            },
        }
        
        return report
    
    # Private analysis methods (placeholders for now)
    def _run_sliding_window_analysis(self, parameters: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """Run sliding window analysis."""
        logger.info("Running sliding window analysis")
        # TODO: Implement sliding window analysis
        return {"mgc_candidates": [], "statistics": {}}
    
    def _run_blast_analysis(self, parameters: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """Run BLAST analysis."""
        logger.info("Running BLAST analysis")
        # TODO: Implement BLAST analysis
        return {"blast_results": [], "statistics": {}}
    
    def _run_phylogenetic_analysis(self, parameters: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """Run phylogenetic analysis."""
        logger.info("Running phylogenetic analysis")
        # TODO: Implement phylogenetic analysis
        return {"tree_data": {}, "statistics": {}}
    
    def _run_ml_analysis(self, parameters: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """Run machine learning analysis."""
        logger.info("Running machine learning analysis")
        # TODO: Implement ML analysis
        return {"predictions": [], "model_metrics": {}}
    
    def _run_pathway_analysis(self, parameters: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """Run pathway analysis."""
        logger.info("Running pathway analysis")
        # TODO: Implement pathway analysis
        return {"pathway_data": {}, "enrichment_results": {}}
    
    def _run_comparative_analysis(self, parameters: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """Run comparative genomics analysis."""
        logger.info("Running comparative genomics analysis")
        # TODO: Implement comparative analysis
        return {"comparative_data": {}, "statistics": {}}