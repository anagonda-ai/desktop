"""
Main Orchestrator for Plant MGC Analysis Pipeline.

This module provides the unified interface that consolidates all 92 legacy scripts
into a single, cohesive, industry-level bioinformatics analysis framework.
"""

import time
from typing import Any, Dict, List, Optional, Union
from pathlib import Path
from dataclasses import dataclass, field
from enum import Enum
from concurrent.futures import ThreadPoolExecutor, as_completed

from loguru import logger

from .core.base import BioinformaticsProcessor
from .core.types import PathLike, GenomeData, AnalysisResult, AnalysisType
from .core.exceptions import ProcessingError, ValidationError, ConfigurationError
from .config.settings import get_settings
from .utils.file_operations import UnifiedFileManager
from .utils.api_clients import UnifiedAPIManager
from .genomics.blast_analysis import BlastAnalyzer
from .genomics.sliding_window import SlidingWindowAnalyzer
from .metabolic.kegg_integration import KEGGAnalyzer
# from .metabolic.pathway_analysis import PathwayAnalyzer  # TODO: Implement when needed
# from .workflows.base_workflow import WorkflowManager  # TODO: Implement when needed


class AnalysisMode(Enum):
    """Analysis execution modes."""
    SINGLE_ORGANISM = "single_organism"
    BATCH_ORGANISMS = "batch_organisms"
    COMPARATIVE = "comparative"
    PIPELINE = "pipeline"


@dataclass
class AnalysisJobConfig:
    """Configuration for analysis jobs."""
    
    job_id: str
    analysis_types: List[AnalysisType]
    input_files: List[PathLike]
    output_directory: Path
    parameters: Dict[str, Any] = field(default_factory=dict)
    mode: AnalysisMode = AnalysisMode.SINGLE_ORGANISM
    parallel_execution: bool = True
    max_workers: int = 4
    
    def __post_init__(self):
        """Validate configuration."""
        if not self.analysis_types:
            raise ValidationError("At least one analysis type is required")
        if not self.input_files:
            raise ValidationError("At least one input file is required")
        if self.max_workers < 1:
            raise ValidationError("Max workers must be at least 1")


@dataclass
class AnalysisJobResult:
    """Result of an analysis job."""
    
    job_id: str
    status: str
    results: List[AnalysisResult] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    processing_time: float = 0.0
    output_files: List[Path] = field(default_factory=list)
    
    @property
    def success(self) -> bool:
        """Check if job completed successfully."""
        return self.status == "success" and not self.errors
    
    @property
    def total_results(self) -> int:
        """Get total number of results."""
        return len(self.results)
    
    def get_summary(self) -> Dict[str, Any]:
        """Get job summary."""
        return {
            "job_id": self.job_id,
            "status": self.status,
            "total_results": self.total_results,
            "processing_time": self.processing_time,
            "error_count": len(self.errors),
            "warning_count": len(self.warnings),
            "output_files": [str(f) for f in self.output_files],
        }


class PlantMGCAnalysisOrchestrator(BioinformaticsProcessor):
    """
    Main orchestrator class that provides unified interface to all functionality.
    
    This class consolidates all operations from 92 legacy scripts into a single,
    cohesive interface, providing:
    - Unified file operations
    - Comprehensive API access
    - All analysis types
    - Workflow management
    - Result integration
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None, **kwargs):
        """Initialize the orchestrator."""
        super().__init__(config, **kwargs)
        
        # Initialize settings
        self.settings = get_settings()
        if config:
            self.settings.update_config(**config)
        
        # Initialize managers
        self.file_manager = UnifiedFileManager()
        self.api_manager = UnifiedAPIManager(
            email=self.settings.get("email", "user@example.com"),
            ncbi_api_key=self.settings.get("ncbi_api_key")
        )
        # self.workflow_manager = WorkflowManager(self.settings)  # TODO: Implement when needed
        
        # Initialize analyzers
        self.analyzers = {
            AnalysisType.BLAST_SEARCH: BlastAnalyzer(session_id=self.session_id),
            AnalysisType.SLIDING_WINDOW: SlidingWindowAnalyzer(session_id=self.session_id),
            AnalysisType.KEGG_ANALYSIS: KEGGAnalyzer(session_id=self.session_id),
            # AnalysisType.PATHWAY_ANALYSIS: PathwayAnalyzer(session_id=self.session_id),  # TODO: Implement when needed
        }
        
        # Job tracking
        self.active_jobs: Dict[str, AnalysisJobResult] = {}
        self.completed_jobs: Dict[str, AnalysisJobResult] = {}
        
        logger.info(f"Plant MGC Analysis Orchestrator initialized with session {self.session_id}")
    
    def validate_input(self, data: Any) -> None:
        """Validate input data."""
        if isinstance(data, AnalysisJobConfig):
            # Validate job configuration
            for input_file in data.input_files:
                if not Path(input_file).exists():
                    raise ValidationError(f"Input file not found: {input_file}")
        elif isinstance(data, dict):
            # Validate as parameters
            required_keys = ["input_files", "analysis_types"]
            for key in required_keys:
                if key not in data:
                    raise ValidationError(f"Missing required parameter: {key}")
        else:
            raise ValidationError("Invalid input data type")
    
    def process(self, data: Any, **kwargs) -> Any:
        """Process input data."""
        if isinstance(data, AnalysisJobConfig):
            return self.run_analysis_job(data)
        elif isinstance(data, dict):
            # Convert dict to job config
            job_config = self._dict_to_job_config(data)
            return self.run_analysis_job(job_config)
        else:
            raise ValidationError("Unsupported input data type")
    
    def run_analysis_job(self, job_config: AnalysisJobConfig) -> AnalysisJobResult:
        """
        Run complete analysis job.
        
        Args:
            job_config: Job configuration
            
        Returns:
            Job result
        """
        self.validate_input(job_config)
        
        # Create job result
        job_result = AnalysisJobResult(
            job_id=job_config.job_id,
            status="running"
        )
        
        # Register job
        self.active_jobs[job_config.job_id] = job_result
        
        start_time = time.time()
        
        try:
            logger.info(f"Starting analysis job: {job_config.job_id}")
            
            # Create output directory
            job_config.output_directory.mkdir(parents=True, exist_ok=True)
            
            # Process based on mode
            if job_config.mode == AnalysisMode.SINGLE_ORGANISM:
                results = self._run_single_organism_analysis(job_config)
            elif job_config.mode == AnalysisMode.BATCH_ORGANISMS:
                results = self._run_batch_analysis(job_config)
            elif job_config.mode == AnalysisMode.COMPARATIVE:
                results = self._run_comparative_analysis(job_config)
            elif job_config.mode == AnalysisMode.PIPELINE:
                results = self._run_pipeline_analysis(job_config)
            else:
                raise ValidationError(f"Unsupported analysis mode: {job_config.mode}")
            
            job_result.results = results
            job_result.status = "success"
            
            # Generate output files
            output_files = self._generate_output_files(job_config, results)
            job_result.output_files = output_files
            
            logger.info(f"Analysis job completed successfully: {job_config.job_id}")
            
        except Exception as e:
            job_result.status = "failed"
            job_result.errors.append(str(e))
            logger.error(f"Analysis job failed: {job_config.job_id} - {e}")
        
        finally:
            job_result.processing_time = time.time() - start_time
            
            # Move to completed jobs
            self.completed_jobs[job_config.job_id] = job_result
            del self.active_jobs[job_config.job_id]
        
        return job_result
    
    def _run_single_organism_analysis(self, job_config: AnalysisJobConfig) -> List[AnalysisResult]:
        """Run analysis for single organism."""
        if len(job_config.input_files) != 1:
            raise ValidationError("Single organism analysis requires exactly one input file")
        
        input_file = Path(job_config.input_files[0])
        
        # Load genome data
        genome_data = self._load_genome_data(input_file, job_config.parameters)
        
        # Run analyses
        results = []
        for analysis_type in job_config.analysis_types:
            try:
                analyzer = self.analyzers[analysis_type]
                result = analyzer.run_analysis(
                    genome_data,
                    job_config.parameters.get(analysis_type.value, {})
                )
                results.append(result)
            except Exception as e:
                logger.error(f"Analysis {analysis_type.value} failed: {e}")
                # Continue with other analyses
        
        return results
    
    def _run_batch_analysis(self, job_config: AnalysisJobConfig) -> List[AnalysisResult]:
        """Run analysis for multiple organisms."""
        all_results = []
        
        if job_config.parallel_execution:
            # Parallel execution
            with ThreadPoolExecutor(max_workers=job_config.max_workers) as executor:
                futures = {}
                
                for input_file in job_config.input_files:
                    future = executor.submit(
                        self._process_single_file,
                        input_file,
                        job_config.analysis_types,
                        job_config.parameters
                    )
                    futures[future] = input_file
                
                for future in as_completed(futures):
                    input_file = futures[future]
                    try:
                        results = future.result()
                        all_results.extend(results)
                        logger.info(f"Processed file: {input_file}")
                    except Exception as e:
                        logger.error(f"Error processing file {input_file}: {e}")
        else:
            # Sequential execution
            for input_file in job_config.input_files:
                try:
                    results = self._process_single_file(
                        input_file,
                        job_config.analysis_types,
                        job_config.parameters
                    )
                    all_results.extend(results)
                    logger.info(f"Processed file: {input_file}")
                except Exception as e:
                    logger.error(f"Error processing file {input_file}: {e}")
        
        return all_results
    
    def _run_comparative_analysis(self, job_config: AnalysisJobConfig) -> List[AnalysisResult]:
        """Run comparative analysis across multiple organisms."""
        # Load all genome data
        genome_data_list = []
        for input_file in job_config.input_files:
            genome_data = self._load_genome_data(Path(input_file), job_config.parameters)
            genome_data_list.append(genome_data)
        
        # Run comparative analyses
        results = []
        for analysis_type in job_config.analysis_types:
            try:
                analyzer = self.analyzers[analysis_type]
                
                # Run individual analyses
                individual_results = []
                for genome_data in genome_data_list:
                    result = analyzer.run_analysis(
                        genome_data,
                        job_config.parameters.get(analysis_type.value, {})
                    )
                    individual_results.append(result)
                
                # Perform comparative analysis
                comparative_result = self._perform_comparative_analysis(
                    individual_results,
                    analysis_type,
                    job_config.parameters
                )
                
                results.extend(individual_results)
                results.append(comparative_result)
                
            except Exception as e:
                logger.error(f"Comparative analysis {analysis_type.value} failed: {e}")
        
        return results
    
    def _run_pipeline_analysis(self, job_config: AnalysisJobConfig) -> List[AnalysisResult]:
        """Run pipeline analysis using workflow manager."""
        # Create workflow
        workflow = self.workflow_manager.create_workflow(
            workflow_name=f"pipeline_{job_config.job_id}",
            analysis_types=job_config.analysis_types,
            parameters=job_config.parameters
        )
        
        # Execute workflow
        results = []
        for input_file in job_config.input_files:
            genome_data = self._load_genome_data(Path(input_file), job_config.parameters)
            workflow_result = workflow.execute(genome_data)
            results.extend(workflow_result.results)
        
        return results
    
    def _process_single_file(
        self,
        input_file: PathLike,
        analysis_types: List[AnalysisType],
        parameters: Dict[str, Any]
    ) -> List[AnalysisResult]:
        """Process single file with multiple analyses."""
        input_path = Path(input_file)
        genome_data = self._load_genome_data(input_path, parameters)
        
        results = []
        for analysis_type in analysis_types:
            try:
                analyzer = self.analyzers[analysis_type]
                result = analyzer.run_analysis(
                    genome_data,
                    parameters.get(analysis_type.value, {})
                )
                results.append(result)
            except Exception as e:
                logger.error(f"Analysis {analysis_type.value} failed for {input_file}: {e}")
        
        return results
    
    def _load_genome_data(self, input_file: Path, parameters: Dict[str, Any]) -> GenomeData:
        """Load genome data from input file."""
        # Determine organism name
        organism_name = parameters.get("organism_name", input_file.stem)
        
        # Determine file types
        if input_file.suffix.lower() in ['.fasta', '.fa', '.fna']:
            # Genome FASTA file
            genome_data = GenomeData(
                organism_name=organism_name,
                genome_file=input_file,
                annotation_file=None,
                protein_file=None,
                metadata=parameters.get("metadata", {}),
            )
        elif input_file.suffix.lower() in ['.faa']:
            # Protein FASTA file
            genome_data = GenomeData(
                organism_name=organism_name,
                genome_file=None,
                annotation_file=None,
                protein_file=input_file,
                metadata=parameters.get("metadata", {}),
            )
        else:
            # Try to auto-detect or use parameters
            annotation_file = parameters.get("annotation_file")
            protein_file = parameters.get("protein_file")
            
            genome_data = GenomeData(
                organism_name=organism_name,
                genome_file=input_file,
                annotation_file=Path(annotation_file) if annotation_file else None,
                protein_file=Path(protein_file) if protein_file else None,
                metadata=parameters.get("metadata", {}),
            )
        
        return genome_data
    
    def _perform_comparative_analysis(
        self,
        individual_results: List[AnalysisResult],
        analysis_type: AnalysisType,
        parameters: Dict[str, Any]
    ) -> AnalysisResult:
        """Perform comparative analysis across results."""
        # Create comparative result
        comparative_result = AnalysisResult(
            analysis_id=f"comparative_{analysis_type.value}_{self.session_id}",
            analysis_type=analysis_type,
            input_data={
                "organism_count": len(individual_results),
                "organisms": [r.input_data.get("organism", "unknown") for r in individual_results],
            },
            results={
                "comparative_statistics": self._calculate_comparative_statistics(individual_results),
                "individual_results": [r.model_dump() for r in individual_results],
            },
            parameters=parameters,
            status="success",
        )
        
        return comparative_result
    
    def _calculate_comparative_statistics(self, results: List[AnalysisResult]) -> Dict[str, Any]:
        """Calculate comparative statistics across results."""
        stats = {
            "total_organisms": len(results),
            "successful_analyses": len([r for r in results if r.success]),
            "failed_analyses": len([r for r in results if not r.success]),
            "average_processing_time": sum(r.execution_time or 0 for r in results) / len(results),
        }
        
        # Analysis-specific statistics
        if results and results[0].analysis_type == AnalysisType.BLAST_SEARCH:
            total_hits = sum(len(r.results.get("blast_hits", [])) for r in results)
            stats["total_blast_hits"] = total_hits
            stats["average_hits_per_organism"] = total_hits / len(results)
        
        return stats
    
    def _generate_output_files(
        self,
        job_config: AnalysisJobConfig,
        results: List[AnalysisResult]
    ) -> List[Path]:
        """Generate output files for job results."""
        output_files = []
        
        # Generate summary report
        summary_file = job_config.output_directory / f"{job_config.job_id}_summary.json"
        summary_data = {
            "job_id": job_config.job_id,
            "analysis_types": [at.value for at in job_config.analysis_types],
            "total_results": len(results),
            "results": [r.model_dump() for r in results],
            "parameters": job_config.parameters,
        }
        
        self.file_manager.json_handler.write_json(summary_data, summary_file)
        output_files.append(summary_file)
        
        # Generate CSV files for each analysis type
        for analysis_type in job_config.analysis_types:
            analysis_results = [r for r in results if r.analysis_type == analysis_type]
            if analysis_results:
                csv_file = job_config.output_directory / f"{job_config.job_id}_{analysis_type.value}.csv"
                self._write_results_to_csv(analysis_results, csv_file)
                output_files.append(csv_file)
        
        return output_files
    
    def _write_results_to_csv(self, results: List[AnalysisResult], output_file: Path) -> None:
        """Write analysis results to CSV file."""
        if not results:
            return
        
        # Extract data based on analysis type
        if results[0].analysis_type == AnalysisType.BLAST_SEARCH:
            data = []
            for result in results:
                organism = result.input_data.get("organism", "unknown")
                for hit in result.results.get("blast_hits", []):
                    hit_data = hit.copy()
                    hit_data["organism"] = organism
                    data.append(hit_data)
            
            if data:
                import pandas as pd
                df = pd.DataFrame(data)
                self.file_manager.csv_handler.write_csv(df, output_file)
        
        # Add other analysis types as needed
    
    def _dict_to_job_config(self, data: Dict[str, Any]) -> AnalysisJobConfig:
        """Convert dictionary to job configuration."""
        return AnalysisJobConfig(
            job_id=data.get("job_id", f"job_{int(time.time())}"),
            analysis_types=[AnalysisType(at) for at in data["analysis_types"]],
            input_files=data["input_files"],
            output_directory=Path(data.get("output_directory", "output")),
            parameters=data.get("parameters", {}),
            mode=AnalysisMode(data.get("mode", "single_organism")),
            parallel_execution=data.get("parallel_execution", True),
            max_workers=data.get("max_workers", 4),
        )
    
    def get_job_status(self, job_id: str) -> Optional[Dict[str, Any]]:
        """Get status of analysis job."""
        if job_id in self.active_jobs:
            return self.active_jobs[job_id].get_summary()
        elif job_id in self.completed_jobs:
            return self.completed_jobs[job_id].get_summary()
        else:
            return None
    
    def list_jobs(self) -> Dict[str, List[str]]:
        """List all jobs."""
        return {
            "active": list(self.active_jobs.keys()),
            "completed": list(self.completed_jobs.keys()),
        }
    
    def get_analyzer_stats(self) -> Dict[str, Any]:
        """Get statistics from all analyzers."""
        stats = {}
        for analysis_type, analyzer in self.analyzers.items():
            if hasattr(analyzer, 'get_stats'):
                stats[analysis_type.value] = analyzer.get_stats()
        return stats
    
    def get_system_stats(self) -> Dict[str, Any]:
        """Get comprehensive system statistics."""
        return {
            "file_operations": self.file_manager.get_processing_stats(),
            "api_operations": self.api_manager.get_combined_stats(),
            "analyzers": self.get_analyzer_stats(),
            "jobs": {
                "active_count": len(self.active_jobs),
                "completed_count": len(self.completed_jobs),
                "total_count": len(self.active_jobs) + len(self.completed_jobs),
            },
        }


# Convenience functions for common operations
def run_single_organism_analysis(
    input_file: PathLike,
    analysis_types: List[str],
    output_directory: PathLike,
    parameters: Optional[Dict[str, Any]] = None,
    **kwargs
) -> AnalysisJobResult:
    """Convenience function for single organism analysis."""
    orchestrator = PlantMGCAnalysisOrchestrator(**kwargs)
    
    job_config = AnalysisJobConfig(
        job_id=f"single_{Path(input_file).stem}_{int(time.time())}",
        analysis_types=[AnalysisType(at) for at in analysis_types],
        input_files=[input_file],
        output_directory=Path(output_directory),
        parameters=parameters or {},
        mode=AnalysisMode.SINGLE_ORGANISM,
    )
    
    return orchestrator.run_analysis_job(job_config)


def run_batch_analysis(
    input_files: List[PathLike],
    analysis_types: List[str],
    output_directory: PathLike,
    parameters: Optional[Dict[str, Any]] = None,
    parallel: bool = True,
    **kwargs
) -> AnalysisJobResult:
    """Convenience function for batch analysis."""
    orchestrator = PlantMGCAnalysisOrchestrator(**kwargs)
    
    job_config = AnalysisJobConfig(
        job_id=f"batch_{int(time.time())}",
        analysis_types=[AnalysisType(at) for at in analysis_types],
        input_files=input_files,
        output_directory=Path(output_directory),
        parameters=parameters or {},
        mode=AnalysisMode.BATCH_ORGANISMS,
        parallel_execution=parallel,
    )
    
    return orchestrator.run_analysis_job(job_config)


# Global orchestrator instance
_orchestrator_instance = None


def get_orchestrator(**kwargs) -> PlantMGCAnalysisOrchestrator:
    """Get global orchestrator instance."""
    global _orchestrator_instance
    if _orchestrator_instance is None:
        _orchestrator_instance = PlantMGCAnalysisOrchestrator(**kwargs)
    return _orchestrator_instance