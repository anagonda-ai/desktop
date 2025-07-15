#!/usr/bin/env python3
"""
State-of-the-art KEGG metabolic pathway extractor for cluster processing.

This module fetches KEGG organism lists, identifies metabolic modules, and
submits SLURM jobs for large-scale organism processing with professional
job management and comprehensive monitoring.

The script:
1. Fetches all plant organisms from KEGG
2. Identifies all metabolic modules
3. Submits SLURM jobs for each organism
4. Monitors job submission and provides comprehensive reporting

Usage:
    python kegg_pathway_extractor_cluster.py [output_folder] [slurm_script_path]

Example:
    python kegg_pathway_extractor_cluster.py /data/metabolic_cluster ./scripts/kegg_download.sh
"""

import sys
import os
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Set, Any, Tuple
from dataclasses import dataclass
import pandas as pd
import logging
import time

# Import our new framework
sys.path.insert(0, str(Path(__file__).parent.parent))
from core import BaseProcessor, Config, get_config
from core.exceptions import ProcessingError, ValidationError, APIError
from core.types import ProcessingConfig, AnalysisResult, ProcessingStatus
from utils import (
    KEGGClient, CSVHandler, FileManager, BatchProcessor, 
    ProgressTracker, RateLimitedExecutor, ExecutorType
)


@dataclass
class ClusterProcessingConfig(ProcessingConfig):
    """Configuration for cluster-based KEGG processing."""
    max_concurrent_jobs: int = 100
    job_submission_delay: float = 1.0  # Seconds between job submissions
    organism_batch_size: int = 50
    module_check_workers: int = 10
    validate_slurm_script: bool = True
    dry_run: bool = False  # For testing without actual job submission
    slurm_partition: str = "general"
    slurm_memory: str = "8G"
    slurm_time: str = "24:00:00"


@dataclass
class OrganismInfo:
    """Information about a KEGG organism."""
    code: str
    name: str
    full_name: str
    is_plant: bool = True
    
    @classmethod
    def from_kegg_line(cls, line: str) -> Optional['OrganismInfo']:
        """Create OrganismInfo from KEGG organism list line."""
        try:
            parts = line.split("\t")
            if len(parts) >= 4:
                org_code = parts[1]
                org_name = parts[2]
                full_name = parts[2] if "(" not in parts[2] else parts[2].split(" (")[0]
                is_plant = "Plants" in parts[3]
                
                return cls(
                    code=org_code,
                    name=org_name,
                    full_name=full_name,
                    is_plant=is_plant
                )
            return None
        except Exception:
            return None


@dataclass
class MetabolicModule:
    """Information about a metabolic module."""
    module_id: str
    name: str
    class_info: str
    is_metabolic: bool = True


@dataclass
class JobSubmissionResult:
    """Result of SLURM job submission."""
    organism_code: str
    job_id: Optional[str] = None
    success: bool = False
    error_message: Optional[str] = None
    command: Optional[str] = None


class SLURMJobManager:
    """
    Professional SLURM job management system.
    
    Handles job submission with rate limiting, monitoring, and
    comprehensive error handling.
    """
    
    def __init__(self, config: ClusterProcessingConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.submitted_jobs: List[JobSubmissionResult] = []
        self.rate_limited_executor = RateLimitedExecutor(
            rate_limit=config.job_submission_delay,
            max_workers=1,  # Sequential job submission
            executor_type=ExecutorType.THREAD
        )
    
    def validate_slurm_script(self, script_path: Path) -> bool:
        """
        Validate SLURM script exists and is executable.
        
        Args:
            script_path: Path to SLURM script
            
        Returns:
            True if script is valid
            
        Raises:
            ValidationError: If script validation fails
        """
        if not script_path.exists():
            raise ValidationError(
                f"SLURM script not found: {script_path}",
                field_name="slurm_script_path"
            )
        
        if not script_path.is_file():
            raise ValidationError(
                f"SLURM script path is not a file: {script_path}",
                field_name="slurm_script_path"
            )
        
        # Check if file is executable
        if not os.access(script_path, os.X_OK):
            self.logger.warning(f"SLURM script is not executable: {script_path}")
            # Try to make it executable
            try:
                script_path.chmod(0o755)
                self.logger.info(f"Made script executable: {script_path}")
            except Exception as e:
                raise ValidationError(
                    f"Cannot make script executable: {e}",
                    field_name="slurm_script_path"
                )
        
        return True
    
    def submit_job(
        self, 
        organism_code: str, 
        script_path: Path, 
        output_folder: Path,
        metabolic_modules_path: Path
    ) -> JobSubmissionResult:
        """
        Submit SLURM job for organism processing.
        
        Args:
            organism_code: KEGG organism code
            script_path: Path to SLURM script
            output_folder: Output folder for results
            metabolic_modules_path: Path to metabolic modules CSV
            
        Returns:
            JobSubmissionResult with submission details
        """
        # Construct sbatch command
        command = (
            f"sbatch "
            f"--partition={self.config.slurm_partition} "
            f"--mem={self.config.slurm_memory} "
            f"--time={self.config.slurm_time} "
            f"--job-name=kegg_{organism_code} "
            f"{script_path} {organism_code} {output_folder} {metabolic_modules_path}"
        )
        
        result = JobSubmissionResult(
            organism_code=organism_code,
            command=command
        )
        
        if self.config.dry_run:
            # Simulate job submission
            result.success = True
            result.job_id = f"dry_run_{organism_code}"
            self.logger.info(f"DRY RUN: Would submit job for {organism_code}")
            return result
        
        try:
            # Submit job
            self.logger.debug(f"Submitting job for {organism_code}")
            process_result = subprocess.run(
                command,
                shell=True,
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if process_result.returncode == 0:
                # Parse job ID from sbatch output
                output = process_result.stdout.strip()
                if "Submitted batch job" in output:
                    job_id = output.split()[-1]
                    result.job_id = job_id
                    result.success = True
                    self.logger.info(f"✅ Submitted job {job_id} for organism {organism_code}")
                else:
                    result.error_message = f"Unexpected sbatch output: {output}"
                    self.logger.error(f"Unexpected sbatch output for {organism_code}: {output}")
            else:
                result.error_message = process_result.stderr.strip()
                self.logger.error(f"Job submission failed for {organism_code}: {result.error_message}")
            
        except subprocess.TimeoutExpired:
            result.error_message = "Job submission timed out"
            self.logger.error(f"Job submission timeout for {organism_code}")
        except Exception as e:
            result.error_message = str(e)
            self.logger.error(f"Job submission error for {organism_code}: {e}")
        
        self.submitted_jobs.append(result)
        return result
    
    def submit_jobs_batch(
        self,
        organisms: List[OrganismInfo],
        script_path: Path,
        output_folder: Path,
        metabolic_modules_path: Path
    ) -> List[JobSubmissionResult]:
        """
        Submit jobs for a batch of organisms with rate limiting.
        
        Args:
            organisms: List of organisms to process
            script_path: Path to SLURM script
            output_folder: Output folder for results
            metabolic_modules_path: Path to metabolic modules CSV
            
        Returns:
            List of job submission results
        """
        results = []
        
        def submit_single_job(organism: OrganismInfo) -> JobSubmissionResult:
            return self.submit_job(organism.code, script_path, output_folder, metabolic_modules_path)
        
        # Submit jobs with rate limiting
        with self.rate_limited_executor as executor:
            job_results = executor.map(submit_single_job, organisms)
            results.extend(job_results)
        
        return results
    
    def get_submission_stats(self) -> Dict[str, Any]:
        """Get comprehensive job submission statistics."""
        total_jobs = len(self.submitted_jobs)
        successful_jobs = sum(1 for job in self.submitted_jobs if job.success)
        failed_jobs = total_jobs - successful_jobs
        
        return {
            "total_submitted": total_jobs,
            "successful": successful_jobs,
            "failed": failed_jobs,
            "success_rate": (successful_jobs / total_jobs * 100) if total_jobs > 0 else 0,
            "job_ids": [job.job_id for job in self.submitted_jobs if job.job_id],
            "failed_organisms": [
                job.organism_code for job in self.submitted_jobs 
                if not job.success
            ]
        }


class KEGGOrganismManager:
    """
    Professional KEGG organism discovery and management.
    
    Fetches and manages KEGG organism information with caching
    and comprehensive validation.
    """
    
    def __init__(self):
        self.kegg_client = KEGGClient()
        self.csv_handler = CSVHandler()
        self.logger = logging.getLogger(__name__)
    
    def fetch_all_organisms(self) -> List[OrganismInfo]:
        """
        Fetch all organisms from KEGG API.
        
        Returns:
            List of OrganismInfo objects
            
        Raises:
            APIError: If KEGG API request fails
        """
        self.logger.info("Fetching KEGG organism list...")
        
        try:
            response = self.kegg_client.get("list/organism")
            
            if response.status_code != 200:
                raise APIError(
                    f"KEGG API returned status {response.status_code}",
                    api_name="KEGG",
                    endpoint="list/organism",
                    status_code=response.status_code
                )
            
            # Parse organisms
            organisms = []
            lines = response.data.strip().split("\n")
            
            for line in lines:
                organism = OrganismInfo.from_kegg_line(line)
                if organism:
                    organisms.append(organism)
            
            self.logger.info(f"Found {len(organisms)} total organisms")
            return organisms
            
        except Exception as e:
            if isinstance(e, APIError):
                raise
            raise APIError(
                f"Failed to fetch organisms: {e}",
                api_name="KEGG",
                endpoint="list/organism",
                original_error=e
            )
    
    def filter_plant_organisms(self, organisms: List[OrganismInfo]) -> List[OrganismInfo]:
        """Filter organisms to include only plants."""
        plants = [org for org in organisms if org.is_plant]
        self.logger.info(f"Found {len(plants)} plant organisms")
        return plants
    
    def save_organisms_csv(self, organisms: List[OrganismInfo], output_path: Path) -> None:
        """
        Save organisms to CSV file.
        
        Args:
            organisms: List of organisms to save
            output_path: Path to output CSV file
        """
        data = [
            {
                "Organism_Code": org.code,
                "Organism_Name": org.name,
                "Full_Name": org.full_name,
                "Is_Plant": org.is_plant
            }
            for org in organisms
        ]
        
        df = pd.DataFrame(data)
        self.csv_handler.write_csv(df, output_path)
        self.logger.info(f"Saved {len(organisms)} organisms to {output_path}")
    
    def load_organisms_csv(self, csv_path: Path) -> List[OrganismInfo]:
        """
        Load organisms from existing CSV file.
        
        Args:
            csv_path: Path to CSV file
            
        Returns:
            List of OrganismInfo objects
        """
        self.logger.info(f"Loading existing organism list from {csv_path}")
        
        df = self.csv_handler.read_csv(csv_path)
        
        organisms = []
        for _, row in df.iterrows():
            organism = OrganismInfo(
                code=row["Organism_Code"],
                name=row["Organism_Name"],
                full_name=row.get("Full_Name", row["Organism_Name"]),
                is_plant=row.get("Is_Plant", True)
            )
            organisms.append(organism)
        
        self.logger.info(f"Loaded {len(organisms)} organisms from CSV")
        return organisms


class MetabolicModuleManager:
    """
    Professional metabolic module discovery and management.
    
    Identifies metabolic modules from KEGG with batch processing
    and comprehensive validation.
    """
    
    def __init__(self, config: ClusterProcessingConfig):
        self.config = config
        self.kegg_client = KEGGClient()
        self.csv_handler = CSVHandler()
        self.logger = logging.getLogger(__name__)
        self.batch_processor = BatchProcessor(
            max_workers=config.module_check_workers,
            executor_type=ExecutorType.THREAD
        )
    
    def fetch_all_modules(self) -> List[str]:
        """
        Fetch all module IDs from KEGG.
        
        Returns:
            List of module IDs
        """
        self.logger.info("Fetching all KEGG modules...")
        
        response = self.kegg_client.get("list/module")
        
        if response.status_code != 200:
            raise APIError(
                f"Failed to fetch modules: status {response.status_code}",
                api_name="KEGG",
                endpoint="list/module"
            )
        
        # Parse module IDs
        module_ids = []
        lines = response.data.strip().splitlines()
        
        for line in lines:
            parts = line.split("\t")
            if parts:
                module_ids.append(parts[0])
        
        self.logger.info(f"Found {len(module_ids)} total modules")
        return module_ids
    
    def check_module_metabolic(self, module_id: str) -> Optional[MetabolicModule]:
        """
        Check if a module is metabolic.
        
        Args:
            module_id: KEGG module ID
            
        Returns:
            MetabolicModule if metabolic, None otherwise
        """
        try:
            response = self.kegg_client.get(f"get/{module_id}")
            
            if response.status_code != 200:
                self.logger.warning(f"Failed to fetch module {module_id}")
                return None
            
            # Parse module information
            lines = response.data.splitlines()
            module_name = ""
            module_class = ""
            
            for line in lines:
                if line.startswith("NAME"):
                    module_name = line[12:].strip()
                elif line.startswith("CLASS"):
                    module_class = line[12:].strip()
            
            # Check if metabolic
            is_metabolic = "metabolism" in module_class.lower()
            
            if is_metabolic:
                return MetabolicModule(
                    module_id=module_id,
                    name=module_name,
                    class_info=module_class,
                    is_metabolic=True
                )
            
            return None
            
        except Exception as e:
            self.logger.error(f"Error checking module {module_id}: {e}")
            return None
    
    def identify_metabolic_modules(self, module_ids: List[str]) -> List[MetabolicModule]:
        """
        Identify metabolic modules from list of module IDs.
        
        Args:
            module_ids: List of module IDs to check
            
        Returns:
            List of metabolic modules
        """
        self.logger.info(f"Checking {len(module_ids)} modules for metabolic classification...")
        
        # Process modules in batches
        def batch_function(module_batch):
            metabolic_modules = []
            for module_id in module_batch:
                module = self.check_module_metabolic(module_id)
                if module:
                    metabolic_modules.append(module)
            return metabolic_modules
        
        batch_results = self.batch_processor.process_parallel(
            module_ids,
            batch_function,
            description="Checking modules for metabolic classification"
        )
        
        # Collect all metabolic modules
        all_metabolic_modules = []
        for batch_result in batch_results:
            if batch_result.success:
                all_metabolic_modules.extend(batch_result.result)
        
        self.logger.info(f"Found {len(all_metabolic_modules)} metabolic modules")
        return all_metabolic_modules
    
    def save_metabolic_modules_csv(self, modules: List[MetabolicModule], output_path: Path) -> None:
        """
        Save metabolic modules to CSV file.
        
        Args:
            modules: List of metabolic modules
            output_path: Path to output CSV file
        """
        data = [
            {
                "Module_ID": module.module_id,
                "Name": module.name,
                "Class": module.class_info
            }
            for module in modules
        ]
        
        df = pd.DataFrame(data)
        self.csv_handler.write_csv(df, output_path)
        self.logger.info(f"Saved {len(modules)} metabolic modules to {output_path}")


class KEGGClusterProcessor(BaseProcessor):
    """
    State-of-the-art KEGG cluster processing orchestrator.
    
    Manages the complete workflow for large-scale KEGG organism processing
    including organism discovery, metabolic module identification, and
    SLURM job submission with comprehensive monitoring.
    """
    
    def __init__(self, config: Optional[ClusterProcessingConfig] = None):
        super().__init__(config)
        self.config = config or ClusterProcessingConfig()
        
        # Initialize managers
        self.file_manager = FileManager()
        self.organism_manager = KEGGOrganismManager()
        self.module_manager = MetabolicModuleManager(self.config)
        self.job_manager = SLURMJobManager(self.config)
    
    def validate_input(self, data: Dict[str, Any]) -> None:
        """
        Validate input parameters.
        
        Args:
            data: Input parameters dictionary
            
        Raises:
            ValidationError: If validation fails
        """
        # Validate output folder
        output_folder = Path(data.get("output_folder", ""))
        if output_folder.exists() and not output_folder.is_dir():
            raise ValidationError(
                f"Output path exists but is not a directory: {output_folder}",
                field_name="output_folder"
            )
        
        # Validate SLURM script
        slurm_script_path = Path(data.get("slurm_script_path", ""))
        if self.config.validate_slurm_script:
            self.job_manager.validate_slurm_script(slurm_script_path)
    
    def setup_output_structure(self, output_folder: Path) -> Dict[str, Path]:
        """
        Setup output directory structure.
        
        Args:
            output_folder: Root output folder
            
        Returns:
            Dictionary of output paths
        """
        paths = {
            "root": output_folder,
            "plants_csv": output_folder / "plants_list.csv",
            "metabolic_modules_csv": output_folder / "metabolic_modules.csv"
        }
        
        # Create directories
        self.file_manager.ensure_directory(output_folder)
        
        return paths
    
    def process_organisms(self, output_paths: Dict[str, Path]) -> List[OrganismInfo]:
        """
        Process organism discovery and management.
        
        Args:
            output_paths: Dictionary of output paths
            
        Returns:
            List of plant organisms
        """
        plants_csv = output_paths["plants_csv"]
        
        if plants_csv.exists():
            # Load existing organisms
            organisms = self.organism_manager.load_organisms_csv(plants_csv)
            plants = self.organism_manager.filter_plant_organisms(organisms)
        else:
            # Fetch organisms from KEGG
            all_organisms = self.organism_manager.fetch_all_organisms()
            plants = self.organism_manager.filter_plant_organisms(all_organisms)
            
            # Save to CSV
            self.organism_manager.save_organisms_csv(plants, plants_csv)
        
        return plants
    
    def process_metabolic_modules(self, output_paths: Dict[str, Path]) -> List[MetabolicModule]:
        """
        Process metabolic module discovery and management.
        
        Args:
            output_paths: Dictionary of output paths
            
        Returns:
            List of metabolic modules
        """
        metabolic_csv = output_paths["metabolic_modules_csv"]
        
        if metabolic_csv.exists():
            # Load existing modules (simplified, just IDs)
            df = self.module_manager.csv_handler.read_csv(metabolic_csv)
            metabolic_modules = [
                MetabolicModule(
                    module_id=row["Module_ID"],
                    name=row.get("Name", ""),
                    class_info=row.get("Class", ""),
                    is_metabolic=True
                )
                for _, row in df.iterrows()
            ]
            self.logger.info(f"Loaded {len(metabolic_modules)} existing metabolic modules")
        else:
            # Discover metabolic modules
            all_module_ids = self.module_manager.fetch_all_modules()
            metabolic_modules = self.module_manager.identify_metabolic_modules(all_module_ids)
            
            # Save to CSV
            self.module_manager.save_metabolic_modules_csv(metabolic_modules, metabolic_csv)
        
        return metabolic_modules
    
    def submit_organism_jobs(
        self,
        plants: List[OrganismInfo],
        slurm_script_path: Path,
        output_folder: Path,
        metabolic_modules_path: Path
    ) -> List[JobSubmissionResult]:
        """
        Submit SLURM jobs for all plant organisms.
        
        Args:
            plants: List of plant organisms
            slurm_script_path: Path to SLURM script
            output_folder: Output folder
            metabolic_modules_path: Path to metabolic modules CSV
            
        Returns:
            List of job submission results
        """
        self.logger.info(f"Submitting SLURM jobs for {len(plants)} plant organisms...")
        
        # Submit jobs in batches with progress tracking
        tracker = ProgressTracker(len(plants), "Submitting SLURM jobs")
        
        results = []
        batch_size = self.config.organism_batch_size
        
        for i in range(0, len(plants), batch_size):
            batch = plants[i:i + batch_size]
            
            batch_results = self.job_manager.submit_jobs_batch(
                batch, slurm_script_path, output_folder, metabolic_modules_path
            )
            
            results.extend(batch_results)
            tracker.update(len(batch))
            
            # Small delay between batches to avoid overwhelming the scheduler
            if i + batch_size < len(plants):
                time.sleep(self.config.job_submission_delay)
        
        return results
    
    def process(self, data: Dict[str, Any], **kwargs) -> Dict[str, Any]:
        """
        Process complete KEGG cluster workflow.
        
        Args:
            data: Input parameters
            **kwargs: Additional processing parameters
            
        Returns:
            Dictionary with processing results
        """
        output_folder = Path(data["output_folder"])
        slurm_script_path = Path(data["slurm_script_path"])
        
        # Setup output structure
        output_paths = self.setup_output_structure(output_folder)
        
        # Process organisms
        plants = self.process_organisms(output_paths)
        
        # Process metabolic modules
        metabolic_modules = self.process_metabolic_modules(output_paths)
        
        # Submit organism processing jobs
        job_results = self.submit_organism_jobs(
            plants,
            slurm_script_path,
            output_folder,
            output_paths["metabolic_modules_csv"]
        )
        
        # Get submission statistics
        job_stats = self.job_manager.get_submission_stats()
        
        # Log summary
        self.logger.info(
            f"✅ Cluster processing setup complete:\n"
            f"   Plant organisms: {len(plants)}\n"
            f"   Metabolic modules: {len(metabolic_modules)}\n"
            f"   Jobs submitted: {job_stats['successful']}/{job_stats['total_submitted']}\n"
            f"   Success rate: {job_stats['success_rate']:.1f}%"
        )
        
        return {
            "status": "completed",
            "output_folder": str(output_folder),
            "plant_organisms_count": len(plants),
            "metabolic_modules_count": len(metabolic_modules),
            "job_submission_stats": job_stats,
            "output_files": {
                "plants_csv": str(output_paths["plants_csv"]),
                "metabolic_modules_csv": str(output_paths["metabolic_modules_csv"])
            }
        }


def main():
    """
    Main entry point for KEGG cluster processor.
    
    Provides command-line interface with professional argument handling
    and comprehensive error reporting.
    """
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger(__name__)
    
    # Parse command line arguments
    if len(sys.argv) == 1:
        # Use default paths for backward compatibility
        output_folder = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic_cluster"
        slurm_script_path = "/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/kegg_organisms_download.sh"
    elif len(sys.argv) == 3:
        output_folder = sys.argv[1]
        slurm_script_path = sys.argv[2]
    else:
        print("Usage: python kegg_pathway_extractor_cluster.py [output_folder] [slurm_script_path]")
        print("\nArguments:")
        print("  output_folder      : Directory for output files and job results")
        print("  slurm_script_path  : Path to SLURM script for organism processing")
        print("\nExample:")
        print("  python kegg_pathway_extractor_cluster.py /data/cluster ./scripts/kegg_download.sh")
        print("\nNote: If no arguments provided, uses default paths for backward compatibility")
        sys.exit(1)
    
    # Create processor with configuration
    config = ClusterProcessingConfig(
        max_concurrent_jobs=100,
        job_submission_delay=1.0,
        organism_batch_size=50,
        module_check_workers=10,
        validate_slurm_script=True,
        dry_run=False  # Set to True for testing
    )
    
    processor = KEGGClusterProcessor(config)
    
    try:
        # Prepare input data
        input_data = {
            "output_folder": output_folder,
            "slurm_script_path": slurm_script_path
        }
        
        # Run processing
        logger.info("Starting KEGG cluster processing workflow")
        result = processor.run(input_data)
        
        if result.status == ProcessingStatus.COMPLETED:
            logger.info("Cluster processing completed successfully")
            if result.results:
                stats = result.results
                print(f"✅ All metabolic module-based KEGG data collection initiated.")
                print(f"   Plant organisms: {stats['plant_organisms_count']}")
                print(f"   Metabolic modules: {stats['metabolic_modules_count']}")
                print(f"   Jobs submitted: {stats['job_submission_stats']['successful']}")
                
                if stats['job_submission_stats']['failed'] > 0:
                    print(f"   Failed submissions: {stats['job_submission_stats']['failed']}")
        else:
            logger.error(f"Cluster processing failed: {result.error_message}")
            sys.exit(1)
    
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()