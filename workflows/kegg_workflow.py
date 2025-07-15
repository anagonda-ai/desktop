"""
KEGG Analysis Workflow orchestrating organism processing, pathway extraction,
and metabolic module identification into a unified pipeline.
"""

import time
from pathlib import Path
from typing import Dict, Any, List, Optional

from .base_workflow import BaseWorkflow
from core.exceptions import ProcessingError, ValidationError
from core.types import ProcessingStatus
from utils import ProgressTracker, FileManager

# Import the processors we created
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from python_scripts.process_organism import OrganismProcessor, OrganismProcessingConfig
from python_scripts.kegg_pathway_extractor_local import KEGGMetabolicExtractor, MetabolicExtractionConfig
from python_scripts.kegg_pathway_extractor_cluster import KEGGClusterProcessor, ClusterProcessingConfig


class KEGGAnalysisWorkflow(BaseWorkflow):
    """
    Unified KEGG analysis workflow orchestrating organism processing,
    pathway extraction, and cluster-based analysis.
    
    This workflow coordinates:
    1. Organism processing for metabolic modules
    2. Local metabolic pathway extraction
    3. Cluster-based processing for large-scale analysis
    4. Result integration and reporting
    """
    
    def __init__(self):
        super().__init__("KEGG")
        
        # Initialize processors
        self.organism_processor = None
        self.metabolic_extractor = None
        self.cluster_processor = None
        
        # Workflow state
        self.organisms_processed = []
        self.metabolic_modules_found = 0
        self.cluster_jobs_submitted = 0
    
    def get_required_parameters(self) -> List[str]:
        """Get required parameters for KEGG workflow."""
        return ["workspace"]
    
    def get_required_subdirectories(self) -> List[str]:
        """Get required subdirectories for KEGG workflow."""
        return [
            "organisms", 
            "metabolic_modules", 
            "cluster_results",
            "integrated_results",
            "reports"
        ]
    
    def validate_parameters(self, params: Dict[str, Any]) -> None:
        """Validate KEGG workflow parameters."""
        super().validate_parameters(params)
        
        # Validate workflow mode
        mode = params.get('mode', 'full')
        if mode not in ['organism', 'extract', 'cluster', 'full']:
            raise ValidationError(
                f"Invalid workflow mode: {mode}. Must be one of: organism, extract, cluster, full"
            )
        
        # Mode-specific validation
        if mode in ['organism', 'full']:
            if 'organism_codes' not in params and 'organisms_file' not in params:
                raise ValidationError(
                    "Organism processing requires either 'organism_codes' list or 'organisms_file'"
                )
        
        if mode in ['extract', 'full']:
            if 'source_folder' not in params:
                # Use default from organisms processing
                output_paths = self.get_output_paths()
                params['source_folder'] = output_paths['organisms']
        
        if mode == 'cluster':
            if 'slurm_script_path' not in params:
                raise ValidationError(
                    "Cluster mode requires 'slurm_script_path' parameter"
                )
    
    def _setup_processors(self, params: Dict[str, Any]) -> None:
        """Initialize processors with appropriate configurations."""
        
        # Organism processor configuration
        organism_config = OrganismProcessingConfig(
            batch_size=params.get('organism_batch_size', 10),
            max_workers=params.get('organism_workers', 4),
            skip_existing=params.get('skip_existing', True)
        )
        self.organism_processor = OrganismProcessor(organism_config)
        
        # Metabolic extractor configuration
        extract_config = MetabolicExtractionConfig(
            max_workers=params.get('extract_workers', 10),
            cache_metabolic_results=params.get('use_cache', True),
            copy_files=params.get('copy_files', True)
        )
        self.metabolic_extractor = KEGGMetabolicExtractor(extract_config)
        
        # Cluster processor configuration
        cluster_config = ClusterProcessingConfig(
            max_concurrent_jobs=params.get('max_cluster_jobs', 100),
            job_submission_delay=params.get('job_delay', 1.0),
            organism_batch_size=params.get('cluster_batch_size', 50),
            dry_run=params.get('dry_run', False)
        )
        self.cluster_processor = KEGGClusterProcessor(cluster_config)
    
    def _load_organisms(self, params: Dict[str, Any]) -> List[str]:
        """Load organism codes from parameters or file."""
        if 'organism_codes' in params:
            return params['organism_codes']
        
        if 'organisms_file' in params:
            organisms_file = Path(params['organisms_file'])
            if not organisms_file.exists():
                raise ValidationError(f"Organisms file not found: {organisms_file}")
            
            # Read organism codes from file
            with open(organisms_file, 'r') as f:
                organisms = [line.strip() for line in f if line.strip()]
            
            self.logger.info(f"Loaded {len(organisms)} organisms from {organisms_file}")
            return organisms
        
        raise ValidationError("No organism source provided")
    
    def _get_metabolic_modules_file(self, params: Dict[str, Any]) -> Path:
        """Get or create metabolic modules file."""
        if 'metabolic_modules_file' in params:
            return Path(params['metabolic_modules_file'])
        
        # Use default location
        output_paths = self.get_output_paths()
        return output_paths['metabolic_modules'] / "metabolic_modules.csv"
    
    def _process_organisms(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Process organisms for metabolic modules."""
        self.logger.info("Starting organism processing phase")
        
        organisms = self._load_organisms(params)
        output_paths = self.get_output_paths()
        metabolic_modules_file = self._get_metabolic_modules_file(params)
        
        # Setup organism output directory
        organisms_dir = output_paths['organisms']
        
        # Track progress
        tracker = ProgressTracker(len(organisms), "Processing organisms")
        
        results = {
            'processed_organisms': [],
            'failed_organisms': [],
            'total_modules': 0,
            'total_genes': 0
        }
        
        for organism_code in organisms:
            try:
                # Process organism
                input_data = {
                    'organism_code': organism_code,
                    'root_folder': str(organisms_dir),
                    'metabolic_modules_path': str(metabolic_modules_file)
                }
                
                result = self.organism_processor.run(input_data)
                
                if result.status == ProcessingStatus.COMPLETED:
                    organism_stats = result.results
                    results['processed_organisms'].append(organism_code)
                    results['total_modules'] += organism_stats.get('modules_processed', 0)
                    results['total_genes'] += organism_stats.get('genes_processed', 0)
                    
                    self.logger.info(
                        f"✅ {organism_code}: "
                        f"{organism_stats.get('modules_processed', 0)} modules, "
                        f"{organism_stats.get('genes_processed', 0)} genes"
                    )
                else:
                    results['failed_organisms'].append(organism_code)
                    self.logger.error(f"Failed to process organism {organism_code}")
                
                tracker.update(1)
                
            except Exception as e:
                results['failed_organisms'].append(organism_code)
                self.logger.error(f"Error processing organism {organism_code}: {e}")
                tracker.mark_failed(1)
        
        self.organisms_processed = results['processed_organisms']
        
        # Log summary
        self.logger.info(
            f"Organism processing complete: "
            f"{len(results['processed_organisms'])} processed, "
            f"{len(results['failed_organisms'])} failed, "
            f"{results['total_modules']} total modules, "
            f"{results['total_genes']} total genes"
        )
        
        return results
    
    def _extract_metabolic_modules(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Extract metabolic modules from processed data."""
        self.logger.info("Starting metabolic module extraction phase")
        
        output_paths = self.get_output_paths()
        
        # Determine source and output folders
        source_folder = Path(params.get('source_folder', output_paths['organisms']))
        output_folder = output_paths['metabolic_modules']
        
        input_data = {
            'source_folder': str(source_folder),
            'output_folder': str(output_folder)
        }
        
        # Run extraction
        result = self.metabolic_extractor.run(input_data)
        
        if result.status == ProcessingStatus.COMPLETED:
            extraction_stats = result.results
            self.metabolic_modules_found = extraction_stats.get('metabolic_modules_found', 0)
            
            self.logger.info(
                f"Metabolic extraction complete: "
                f"{extraction_stats.get('files_processed', 0)} files processed, "
                f"{self.metabolic_modules_found} metabolic modules found, "
                f"{extraction_stats.get('files_copied', 0)} files copied"
            )
            
            return extraction_stats
        else:
            raise ProcessingError(f"Metabolic extraction failed: {result.error_message}")
    
    def _run_cluster_processing(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Run cluster-based processing for large-scale analysis."""
        self.logger.info("Starting cluster processing phase")
        
        output_paths = self.get_output_paths()
        cluster_output = output_paths['cluster_results']
        
        # Required parameters for cluster processing
        required_cluster_params = ['slurm_script_path']
        for param in required_cluster_params:
            if param not in params:
                raise ValidationError(f"Cluster processing requires '{param}' parameter")
        
        input_data = {
            'output_folder': str(cluster_output),
            'slurm_script_path': params['slurm_script_path']
        }
        
        # Run cluster processing
        result = self.cluster_processor.run(input_data)
        
        if result.status == ProcessingStatus.COMPLETED:
            cluster_stats = result.results
            self.cluster_jobs_submitted = cluster_stats.get('job_submission_stats', {}).get('successful', 0)
            
            self.logger.info(
                f"Cluster processing setup complete: "
                f"{cluster_stats.get('plant_organisms_count', 0)} organisms, "
                f"{cluster_stats.get('metabolic_modules_count', 0)} metabolic modules, "
                f"{self.cluster_jobs_submitted} jobs submitted"
            )
            
            return cluster_stats
        else:
            raise ProcessingError(f"Cluster processing failed: {result.error_message}")
    
    def _integrate_results(self, phase_results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """Integrate results from all phases into unified report."""
        self.logger.info("Integrating workflow results")
        
        integrated_results = {
            'workflow_summary': {
                'total_organisms_processed': len(self.organisms_processed),
                'metabolic_modules_found': self.metabolic_modules_found,
                'cluster_jobs_submitted': self.cluster_jobs_submitted
            },
            'phase_results': phase_results,
            'output_paths': {k: str(v) for k, v in self.get_output_paths().items()}
        }
        
        # Save integrated results
        output_paths = self.get_output_paths()
        results_file = output_paths['integrated_results'] / "kegg_workflow_results.json"
        
        import json
        with open(results_file, 'w') as f:
            json.dump(integrated_results, f, indent=2, default=str)
        
        self.logger.info(f"Saved integrated results to {results_file}")
        
        return integrated_results
    
    def _execute_workflow(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Execute the KEGG analysis workflow."""
        mode = params.get('mode', 'full')
        
        # Setup processors
        self._setup_processors(params)
        
        # Initialize phase results
        phase_results = {}
        
        # Execute workflow phases based on mode
        if mode in ['organism', 'full']:
            phase_results['organism_processing'] = self._process_organisms(params)
        
        if mode in ['extract', 'full']:
            phase_results['metabolic_extraction'] = self._extract_metabolic_modules(params)
        
        if mode in ['cluster', 'full']:
            phase_results['cluster_processing'] = self._run_cluster_processing(params)
        
        # Integrate results
        integrated_results = self._integrate_results(phase_results)
        
        return integrated_results
    
    def pre_execution_hook(self, params: Dict[str, Any]) -> None:
        """Pre-execution setup for KEGG workflow."""
        self.logger.info(f"KEGG workflow mode: {params.get('mode', 'full')}")
        
        # Validate metabolic modules file exists or can be created
        metabolic_modules_file = self._get_metabolic_modules_file(params)
        if not metabolic_modules_file.exists():
            # Create default metabolic modules file if needed
            output_paths = self.get_output_paths()
            
            # Check if we have a default modules file in the project
            default_modules = self.workspace / "data" / "metabolic_modules.csv"
            if default_modules.exists():
                import shutil
                shutil.copy2(default_modules, metabolic_modules_file)
                self.logger.info(f"Copied default metabolic modules to {metabolic_modules_file}")
            else:
                self.logger.warning(f"Metabolic modules file not found: {metabolic_modules_file}")
    
    def post_execution_hook(self, result) -> None:
        """Post-execution cleanup and reporting for KEGG workflow."""
        # Generate summary report
        self._generate_summary_report(result)
    
    def _generate_summary_report(self, result) -> None:
        """Generate a comprehensive summary report."""
        output_paths = self.get_output_paths()
        report_file = output_paths['reports'] / "kegg_workflow_summary.md"
        
        with open(report_file, 'w') as f:
            f.write("# KEGG Analysis Workflow Summary\n\n")
            f.write(f"**Execution Time:** {result.execution_time:.2f} seconds\n")
            f.write(f"**Status:** {result.status.value}\n")
            f.write(f"**Artifacts Generated:** {len(result.artifacts)}\n\n")
            
            if result.results:
                summary = result.results.get('workflow_summary', {})
                f.write("## Results Summary\n\n")
                f.write(f"- **Organisms Processed:** {summary.get('total_organisms_processed', 0)}\n")
                f.write(f"- **Metabolic Modules Found:** {summary.get('metabolic_modules_found', 0)}\n")
                f.write(f"- **Cluster Jobs Submitted:** {summary.get('cluster_jobs_submitted', 0)}\n\n")
            
            # Add metrics
            if result.metrics:
                f.write("## Performance Metrics\n\n")
                for key, value in result.metrics.items():
                    f.write(f"- **{key.replace('_', ' ').title()}:** {value}\n")
                f.write("\n")
            
            # Add output paths
            f.write("## Output Locations\n\n")
            for path_name, path in self.get_output_paths().items():
                f.write(f"- **{path_name.replace('_', ' ').title()}:** `{path}`\n")
        
        self.logger.info(f"Generated summary report: {report_file}")