"""
Industry-Level MGC Comparison Pipeline with SLURM Job Management.

This module provides object-oriented MGC (Metabolic Gene Cluster) comparison pipeline
with comprehensive SLURM job management, error handling, and parallel processing.
"""

import os
import subprocess
import time
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from enum import Enum
import json

from ..core.base_classes import (
    BaseProcessor,
    BioinformaticsConfig,
    ProcessingResult,
    ParallelProcessor
)


class JobStatus(Enum):
    """SLURM job status enumeration."""
    
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"
    TIMEOUT = "TIMEOUT"
    UNKNOWN = "UNKNOWN"


@dataclass
class SLURMJobInfo:
    """Information about a SLURM job."""
    
    job_id: str
    job_name: str
    status: JobStatus = JobStatus.UNKNOWN
    submit_time: Optional[float] = None
    start_time: Optional[float] = None
    end_time: Optional[float] = None
    output_file: Optional[Path] = None
    error_file: Optional[Path] = None
    expected_output: Optional[Path] = None
    
    @property
    def runtime(self) -> Optional[float]:
        """Get job runtime in seconds."""
        if self.start_time and self.end_time:
            return self.end_time - self.start_time
        return None
    
    @property
    def is_finished(self) -> bool:
        """Check if job is finished."""
        return self.status in [JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELLED, JobStatus.TIMEOUT]


@dataclass
class MGCProcessingStep:
    """Configuration for an MGC processing step."""
    
    name: str
    script_path: Path
    job_name_prefix: str
    output_suffix: str
    error_suffix: str
    depends_on: Optional[str] = None
    timeout: int = 3600  # 1 hour default
    
    def get_job_name(self, mgc_name: str) -> str:
        """Get full job name for MGC."""
        return f"{self.job_name_prefix}_{mgc_name}_{self.name}"
    
    def get_output_files(self, mgc_name: str, jobs_dir: Path) -> Tuple[Path, Path]:
        """Get output and error file paths."""
        output_file = jobs_dir / f"{mgc_name}_{self.output_suffix}"
        error_file = jobs_dir / f"{mgc_name}_{self.error_suffix}"
        return output_file, error_file


@dataclass
class PipelineConfig:
    """Configuration for MGC comparison pipeline."""
    
    root_dir: Path
    jobs_dir: Path
    max_concurrent_jobs: int = 100
    max_parallel_mgcs: int = 30
    job_check_interval: int = 30
    file_check_interval: int = 5
    default_timeout: int = 3600
    slurm_user: str = "alongonda"
    skip_existing: bool = True


class SLURMJobManager:
    """
    Professional SLURM job management with monitoring and error handling.
    """
    
    def __init__(self, 
                 config: PipelineConfig,
                 logger):
        """
        Initialize SLURM job manager.
        
        Args:
            config: Pipeline configuration
            logger: Logger instance
        """
        self.config = config
        self.logger = logger
        self.active_jobs: Dict[str, SLURMJobInfo] = {}
        
        # Statistics
        self.job_stats = {
            "jobs_submitted": 0,
            "jobs_completed": 0,
            "jobs_failed": 0,
            "total_wait_time": 0.0,
        }
    
    def wait_for_job_slot(self) -> None:
        """Wait for available job slot based on running jobs."""
        while True:
            running_jobs = self._count_running_jobs()
            
            if running_jobs < self.config.max_concurrent_jobs:
                self.logger.debug(f"Job slot available", running_jobs=running_jobs)
                break
            
            self.logger.info(f"Waiting for job slot", 
                           running_jobs=running_jobs,
                           max_jobs=self.config.max_concurrent_jobs)
            
            time.sleep(self.config.job_check_interval)
    
    def submit_job(self, 
                   job_info: SLURMJobInfo,
                   script_path: Path,
                   mgc_path: Path) -> Optional[str]:
        """
        Submit SLURM job.
        
        Args:
            job_info: Job information
            script_path: Path to SLURM script
            mgc_path: Path to MGC directory
            
        Returns:
            Job ID if successful, None otherwise
        """
        try:
            # Build sbatch command
            cmd = [
                "sbatch",
                "--job-name", job_info.job_name,
                "--output", str(job_info.output_file),
                "--error", str(job_info.error_file),
                str(script_path),
                str(mgc_path)
            ]
            
            # Submit job
            result = subprocess.run(cmd, 
                                  capture_output=True, 
                                  text=True,
                                  timeout=30)
            
            if result.returncode != 0:
                self.logger.error(f"Failed to submit job", 
                                job_name=job_info.job_name,
                                error=result.stderr)
                return None
            
            # Extract job ID
            job_id = result.stdout.strip().split()[-1]
            job_info.job_id = job_id
            job_info.submit_time = time.time()
            job_info.status = JobStatus.PENDING
            
            # Track job
            self.active_jobs[job_id] = job_info
            self.job_stats["jobs_submitted"] += 1
            
            self.logger.info(f"Job submitted successfully", 
                           job_name=job_info.job_name,
                           job_id=job_id)
            
            return job_id
            
        except subprocess.TimeoutExpired:
            self.logger.error(f"Job submission timeout", job_name=job_info.job_name)
            return None
        except Exception as e:
            self.logger.error(f"Job submission failed", 
                            job_name=job_info.job_name,
                            error=str(e))
            return None
    
    def wait_for_job_completion(self, job_id: str, timeout: Optional[int] = None) -> JobStatus:
        """
        Wait for job completion.
        
        Args:
            job_id: SLURM job ID
            timeout: Maximum wait time in seconds
            
        Returns:
            Final job status
        """
        if job_id not in self.active_jobs:
            self.logger.warning(f"Job not tracked", job_id=job_id)
            return JobStatus.UNKNOWN
        
        job_info = self.active_jobs[job_id]
        start_wait_time = time.time()
        timeout = timeout or self.config.default_timeout
        
        self.logger.info(f"Waiting for job completion", 
                       job_id=job_id,
                       job_name=job_info.job_name)
        
        while time.time() - start_wait_time < timeout:
            status = self._check_job_status(job_id)
            job_info.status = status
            
            if status.is_finished:
                wait_time = time.time() - start_wait_time
                self.job_stats["total_wait_time"] += wait_time
                
                if status == JobStatus.COMPLETED:
                    self.job_stats["jobs_completed"] += 1
                    job_info.end_time = time.time()
                    
                    # Wait for output file
                    if job_info.expected_output:
                        self._wait_for_output_file(job_info.expected_output)
                    
                    self.logger.info(f"Job completed successfully", 
                                   job_id=job_id,
                                   wait_time=wait_time)
                else:
                    self.job_stats["jobs_failed"] += 1
                    self.logger.error(f"Job failed", 
                                    job_id=job_id,
                                    status=status.value)
                
                return status
            
            time.sleep(self.config.job_check_interval)
        
        # Timeout reached
        job_info.status = JobStatus.TIMEOUT
        self.job_stats["jobs_failed"] += 1
        
        self.logger.error(f"Job wait timeout", 
                        job_id=job_id,
                        timeout=timeout)
        
        return JobStatus.TIMEOUT
    
    def _count_running_jobs(self) -> int:
        """Count currently running jobs for user."""
        try:
            result = subprocess.run(
                ["squeue", "-u", self.config.slurm_user],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode != 0:
                return 0
            
            lines = result.stdout.strip().splitlines()
            mgc_lines = [line for line in lines if "mgc_" in line]
            
            # Subtract header line if present
            count = len(mgc_lines) - 1 if len(mgc_lines) > 1 else len(mgc_lines)
            return max(0, count)
            
        except Exception as e:
            self.logger.warning(f"Failed to count running jobs", error=str(e))
            return 0
    
    def _check_job_status(self, job_id: str) -> JobStatus:
        """Check status of specific job."""
        try:
            result = subprocess.run(
                ["squeue", "-j", job_id],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode != 0:
                return JobStatus.UNKNOWN
            
            if job_id not in result.stdout:
                # Job not in queue, assume completed
                return JobStatus.COMPLETED
            
            # Parse status from output
            lines = result.stdout.strip().splitlines()
            for line in lines[1:]:  # Skip header
                parts = line.split()
                if len(parts) >= 5 and parts[0] == job_id:
                    status_str = parts[4]
                    
                    if status_str in ["PD"]:
                        return JobStatus.PENDING
                    elif status_str in ["R", "CG"]:
                        return JobStatus.RUNNING
                    elif status_str in ["CD"]:
                        return JobStatus.COMPLETED
                    elif status_str in ["F", "NF"]:
                        return JobStatus.FAILED
                    elif status_str in ["CA"]:
                        return JobStatus.CANCELLED
                    elif status_str in ["TO"]:
                        return JobStatus.TIMEOUT
            
            return JobStatus.UNKNOWN
            
        except Exception as e:
            self.logger.warning(f"Failed to check job status", 
                              job_id=job_id,
                              error=str(e))
            return JobStatus.UNKNOWN
    
    def _wait_for_output_file(self, output_file: Path) -> None:
        """Wait for output file to be ready."""
        self.logger.debug(f"Waiting for output file", path=str(output_file))
        
        # Wait for file existence
        while not output_file.exists():
            time.sleep(self.config.file_check_interval)
        
        # Wait for file to stabilize (no size changes)
        prev_size = -1
        stable_count = 0
        
        while stable_count < 3:  # Require 3 consecutive stable checks
            try:
                curr_size = output_file.stat().st_size
                
                if curr_size == prev_size:
                    stable_count += 1
                else:
                    stable_count = 0
                
                prev_size = curr_size
                time.sleep(self.config.file_check_interval)
                
            except OSError:
                # File might be temporarily unavailable
                stable_count = 0
                time.sleep(self.config.file_check_interval)
        
        self.logger.debug(f"Output file ready", path=str(output_file))
    
    def get_job_statistics(self) -> Dict[str, Any]:
        """Get job management statistics."""
        active_count = len([j for j in self.active_jobs.values() if not j.is_finished])
        
        return {
            **self.job_stats,
            "active_jobs": active_count,
            "tracked_jobs": len(self.active_jobs),
            "average_wait_time": (
                self.job_stats["total_wait_time"] / max(1, self.job_stats["jobs_completed"])
            ),
        }


class MGCComparisonPipeline(ParallelProcessor):
    """
    Professional MGC comparison pipeline with SLURM job orchestration.
    """
    
    def __init__(self, 
                 config: BioinformaticsConfig,
                 pipeline_config: PipelineConfig):
        """
        Initialize MGC comparison pipeline.
        
        Args:
            config: Base processing configuration
            pipeline_config: Pipeline-specific configuration
        """
        super().__init__(config, name="MGCComparisonPipeline")
        
        self.pipeline_config = pipeline_config
        self.job_manager = SLURMJobManager(pipeline_config, self.logger)
        
        # Define processing steps
        self.processing_steps = [
            MGCProcessingStep(
                name="homolog",
                script_path=Path("/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/submit_find_homolog_genes_in_dir.sh"),
                job_name_prefix="mgc",
                output_suffix="homolog.out",
                error_suffix="homolog.err"
            ),
            MGCProcessingStep(
                name="cluster",
                script_path=Path("/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/submit_find_clusters_in_chromosome.sh"),
                job_name_prefix="mgc",
                output_suffix="cluster.out",
                error_suffix="cluster.err",
                depends_on="homolog"
            ),
            MGCProcessingStep(
                name="multi",
                script_path=Path("/groups/itay_mayrose/alongonda/desktop/sh_scripts/powerslurm/daily_usage/submit_multichromosome_statistics.sh"),
                job_name_prefix="mgc",
                output_suffix="multi.out",
                error_suffix="multi.err",
                depends_on="cluster"
            ),
        ]
        
        # Statistics
        self.pipeline_stats = {
            "mgcs_found": 0,
            "mgcs_processed": 0,
            "mgcs_skipped": 0,
            "total_steps_completed": 0,
            "total_processing_time": 0.0,
        }
    
    def process(self, input_data: Any = None, **kwargs) -> ProcessingResult:
        """
        Process all MGC directories in the pipeline.
        
        Args:
            input_data: Optional input data (not used)
            **kwargs: Additional processing parameters
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            self.logger.info(f"Starting MGC comparison pipeline",
                           root_dir=str(self.pipeline_config.root_dir))
            
            # Find MGC directories
            mgc_dirs = self._find_mgc_directories()
            self.pipeline_stats["mgcs_found"] = len(mgc_dirs)
            
            if not mgc_dirs:
                result.add_error("No MGC directories found")
                return result
            
            self.logger.info(f"Found MGC directories", count=len(mgc_dirs))
            
            # Process MGCs in parallel
            mgc_results = self._process_mgcs_parallel(mgc_dirs)
            
            # Aggregate results
            successful_mgcs = []
            total_steps = 0
            
            for mgc_result in mgc_results:
                if mgc_result.success:
                    successful_mgcs.append(mgc_result.metadata.get("mgc_name"))
                    total_steps += mgc_result.metadata.get("steps_completed", 0)
                    
                    if mgc_result.metadata.get("status") == "skipped":
                        self.pipeline_stats["mgcs_skipped"] += 1
                    else:
                        self.pipeline_stats["mgcs_processed"] += 1
                else:
                    result.errors.extend(mgc_result.errors)
                    result.warnings.extend(mgc_result.warnings)
            
            # Update statistics
            self.pipeline_stats["total_steps_completed"] = total_steps
            
            result.metadata.update({
                "mgcs_found": len(mgc_dirs),
                "mgcs_processed": len(successful_mgcs),
                "mgcs_skipped": self.pipeline_stats["mgcs_skipped"],
                "successful_mgcs": successful_mgcs,
                "total_steps_completed": total_steps,
                "pipeline_stats": self.pipeline_stats,
                "job_stats": self.job_manager.get_job_statistics()
            })
            
            if not successful_mgcs:
                result.add_error("No MGCs were processed successfully")
            
            self.logger.info(f"MGC comparison pipeline completed",
                           processed=len(successful_mgcs),
                           total_steps=total_steps)
            
        except Exception as e:
            result.add_error(f"Pipeline failed: {str(e)}")
        
        return result
    
    def _find_mgc_directories(self) -> List[Path]:
        """Find all MGC subdirectories."""
        try:
            mgc_dirs = []
            
            for item in self.pipeline_config.root_dir.iterdir():
                if item.is_dir():
                    mgc_dirs.append(item)
            
            self.logger.debug(f"Found MGC directories", count=len(mgc_dirs))
            return mgc_dirs
            
        except Exception as e:
            self.logger.error(f"Failed to find MGC directories", error=str(e))
            raise
    
    def _process_mgcs_parallel(self, mgc_dirs: List[Path]) -> List[ProcessingResult]:
        """
        Process MGCs in parallel.
        
        Args:
            mgc_dirs: List of MGC directories
            
        Returns:
            List of processing results
        """
        if self.pipeline_config.max_parallel_mgcs == 1:
            # Sequential processing
            return [self._process_single_mgc(mgc_dir) for mgc_dir in mgc_dirs]
        
        # Parallel processing
        results = []
        
        with ThreadPoolExecutor(max_workers=self.pipeline_config.max_parallel_mgcs) as executor:
            future_to_mgc = {
                executor.submit(self._process_single_mgc, mgc_dir): mgc_dir
                for mgc_dir in mgc_dirs
            }
            
            for future in as_completed(future_to_mgc):
                mgc_dir = future_to_mgc[future]
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    error_result = ProcessingResult(success=False)
                    error_result.add_error(f"MGC processing failed for {mgc_dir}: {str(e)}")
                    results.append(error_result)
        
        return results
    
    def _process_single_mgc(self, mgc_dir: Path) -> ProcessingResult:
        """
        Process single MGC through all pipeline steps.
        
        Args:
            mgc_dir: MGC directory path
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        start_time = time.time()
        
        try:
            mgc_name = mgc_dir.name
            
            self.logger.info(f"Processing MGC", mgc_name=mgc_name)
            
            # Ensure jobs directory exists
            self.pipeline_config.jobs_dir.mkdir(parents=True, exist_ok=True)
            
            steps_completed = 0
            
            # Process each step sequentially
            for step in self.processing_steps:
                step_result = self._process_mgc_step(mgc_dir, step)
                
                if not step_result.success:
                    result.errors.extend(step_result.errors)
                    result.warnings.extend(step_result.warnings)
                    break
                
                if step_result.metadata.get("status") == "skipped":
                    self.logger.debug(f"Step skipped", 
                                    mgc_name=mgc_name,
                                    step=step.name)
                else:
                    steps_completed += 1
                    self.logger.info(f"Step completed", 
                                   mgc_name=mgc_name,
                                   step=step.name)
            
            processing_time = time.time() - start_time
            self.pipeline_stats["total_processing_time"] += processing_time
            
            result.metadata.update({
                "mgc_name": mgc_name,
                "mgc_path": str(mgc_dir),
                "steps_completed": steps_completed,
                "total_steps": len(self.processing_steps),
                "processing_time": processing_time,
                "status": "completed" if steps_completed == len(self.processing_steps) else "partial"
            })
            
            if steps_completed == len(self.processing_steps):
                self.logger.info(f"MGC processing completed", 
                               mgc_name=mgc_name,
                               steps=steps_completed)
            else:
                result.add_warning(f"Only {steps_completed}/{len(self.processing_steps)} steps completed")
            
        except Exception as e:
            result.add_error(f"MGC processing failed: {str(e)}")
        
        return result
    
    def _process_mgc_step(self, mgc_dir: Path, step: MGCProcessingStep) -> ProcessingResult:
        """
        Process single MGC step.
        
        Args:
            mgc_dir: MGC directory path
            step: Processing step configuration
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            mgc_name = mgc_dir.name
            output_file, error_file = step.get_output_files(mgc_name, self.pipeline_config.jobs_dir)
            
            # Check if step already completed
            if self.pipeline_config.skip_existing and output_file.exists():
                self.logger.debug(f"Step output exists, skipping", 
                                mgc_name=mgc_name,
                                step=step.name)
                result.metadata["status"] = "skipped"
                return result
            
            # Wait for job slot
            self.job_manager.wait_for_job_slot()
            
            # Create job info
            job_info = SLURMJobInfo(
                job_id="",
                job_name=step.get_job_name(mgc_name),
                output_file=output_file,
                error_file=error_file,
                expected_output=output_file
            )
            
            # Submit job
            job_id = self.job_manager.submit_job(job_info, step.script_path, mgc_dir)
            
            if not job_id:
                result.add_error(f"Failed to submit job for step {step.name}")
                return result
            
            # Wait for completion
            final_status = self.job_manager.wait_for_job_completion(job_id, step.timeout)
            
            if final_status != JobStatus.COMPLETED:
                result.add_error(f"Job failed with status: {final_status.value}")
                return result
            
            result.metadata.update({
                "step": step.name,
                "job_id": job_id,
                "output_file": str(output_file),
                "status": "completed"
            })
            
        except Exception as e:
            result.add_error(f"Step processing failed: {str(e)}")
        
        return result
    
    def get_comprehensive_statistics(self) -> Dict[str, Any]:
        """Get comprehensive pipeline statistics."""
        base_stats = self.get_statistics()
        job_stats = self.job_manager.get_job_statistics()
        
        return {
            **base_stats,
            **self.pipeline_stats,
            "job_management_stats": job_stats,
            "pipeline_config": {
                "root_dir": str(self.pipeline_config.root_dir),
                "jobs_dir": str(self.pipeline_config.jobs_dir),
                "max_concurrent_jobs": self.pipeline_config.max_concurrent_jobs,
                "max_parallel_mgcs": self.pipeline_config.max_parallel_mgcs,
            }
        }
    
    def validate_input(self, input_data: Any = None) -> bool:
        """
        Validate pipeline configuration.
        
        Args:
            input_data: Input data (not used)
            
        Returns:
            True if configuration is valid
        """
        return (self.pipeline_config.root_dir.exists() and
                self.pipeline_config.root_dir.is_dir())


def main():
    """Main entry point for MGC comparison pipeline."""
    # Configuration
    root_dir = Path("/groups/itay_mayrose/alongonda/Plant_MGC/test/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mgc_candidates_dir_fixed")
    jobs_dir = Path("/groups/itay_mayrose/alongonda/desktop/example_jobs")
    
    # Setup configuration
    config = BioinformaticsConfig(
        input_dir=root_dir,
        output_dir=jobs_dir,
        max_workers=1,  # Sequential MGC processing within each parallel thread
        log_level="INFO"
    )
    
    pipeline_config = PipelineConfig(
        root_dir=root_dir,
        jobs_dir=jobs_dir,
        max_concurrent_jobs=100,
        max_parallel_mgcs=30,
        job_check_interval=30,
        skip_existing=True
    )
    
    # Initialize pipeline
    try:
        pipeline = MGCComparisonPipeline(
            config=config,
            pipeline_config=pipeline_config
        )
        
        # Run pipeline
        result = pipeline.run()
        
        # Display results
        if result.success:
            print("🎉 MGC comparison pipeline completed successfully!")
            print(f"   MGCs found: {result.metadata.get('mgcs_found', 0)}")
            print(f"   MGCs processed: {result.metadata.get('mgcs_processed', 0)}")
            print(f"   MGCs skipped: {result.metadata.get('mgcs_skipped', 0)}")
            print(f"   Total steps completed: {result.metadata.get('total_steps_completed', 0)}")
            print(f"   Processing time: {result.processing_time:.2f} seconds")
        else:
            print("❌ MGC comparison pipeline failed!")
            for error in result.errors:
                print(f"   Error: {error}")
        
        # Display warnings
        for warning in result.warnings:
            print(f"   Warning: {warning}")
        
        # Display statistics
        stats = pipeline.get_comprehensive_statistics()
        print(f"\nPipeline Statistics:")
        print(f"   Jobs submitted: {stats['job_management_stats']['jobs_submitted']}")
        print(f"   Jobs completed: {stats['job_management_stats']['jobs_completed']}")
        print(f"   Jobs failed: {stats['job_management_stats']['jobs_failed']}")
        print(f"   Average wait time: {stats['job_management_stats']['average_wait_time']:.2f}s")
        
    except Exception as e:
        print(f"❌ Critical error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())