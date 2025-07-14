"""
Industry-Level KEGG Pathway Extractor with Metabolic Module Detection.

This module provides object-oriented KEGG pathway extraction with comprehensive
error handling, parallel processing, caching, and performance optimization.
"""

import os
import shutil
from typing import Dict, List, Optional, Set, Any, Tuple
from dataclasses import dataclass
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

import requests

from core.base_classes import (
    BaseProcessor,
    BioinformaticsConfig,
    ProcessingResult,
    FileSystemProcessor
)


@dataclass
class MetabolicModuleInfo:
    """Information about a metabolic module."""
    
    module_id: str
    name: str = ""
    classification: str = ""
    is_metabolic: bool = False
    file_path: Optional[Path] = None


@dataclass
class ExtractionConfig:
    """Configuration for KEGG pathway extraction."""
    
    source_folder: Path
    target_folder: Path
    rate_limit_delay: float = 0.35
    max_retries: int = 3
    timeout: int = 30
    copy_files: bool = True
    verify_metabolic: bool = True


class KEGGModuleAPIClient:
    """
    API client for KEGG module information retrieval.
    """
    
    def __init__(self, 
                 base_url: str = "https://rest.kegg.jp",
                 rate_limit_delay: float = 0.35,
                 max_retries: int = 3,
                 timeout: int = 30):
        """
        Initialize KEGG module API client.
        
        Args:
            base_url: KEGG API base URL
            rate_limit_delay: Delay between requests (seconds)
            max_retries: Maximum retry attempts
            timeout: Request timeout (seconds)
        """
        self.base_url = base_url
        self.rate_limit_delay = rate_limit_delay
        self.max_retries = max_retries
        self.timeout = timeout
        
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'KEGGPathwayExtractor/1.0'
        })
        
        self._last_request_time = 0.0
        self._module_cache: Dict[str, bool] = {}
    
    def _enforce_rate_limit(self) -> None:
        """Enforce rate limiting between requests."""
        current_time = time.time()
        time_since_last = current_time - self._last_request_time
        
        if time_since_last < self.rate_limit_delay:
            sleep_time = self.rate_limit_delay - time_since_last
            time.sleep(sleep_time)
        
        self._last_request_time = time.time()
    
    def is_metabolic_module(self, module_id: str) -> bool:
        """
        Check if module is metabolic by querying KEGG API.
        
        Args:
            module_id: KEGG module ID
            
        Returns:
            True if module is metabolic
        """
        # Check cache first
        if module_id in self._module_cache:
            return self._module_cache[module_id]
        
        for attempt in range(self.max_retries):
            try:
                self._enforce_rate_limit()
                
                url = f"{self.base_url}/get/{module_id}"
                response = self.session.get(url, timeout=self.timeout)
                response.raise_for_status()
                
                # Check if module classification contains "metabol"
                is_metabolic = any(
                    "metabol" in line.lower() 
                    for line in response.text.splitlines() 
                    if line.startswith("CLASS")
                )
                
                # Cache result
                self._module_cache[module_id] = is_metabolic
                return is_metabolic
                
            except requests.RequestException as e:
                if attempt == self.max_retries - 1:
                    # Cache as non-metabolic on final failure
                    self._module_cache[module_id] = False
                    return False
                
                # Exponential backoff
                time.sleep(2 ** attempt)
        
        return False
    
    def get_cache_stats(self) -> Dict[str, Any]:
        """Get cache statistics."""
        return {
            "cache_size": len(self._module_cache),
            "metabolic_modules": sum(self._module_cache.values()),
            "non_metabolic_modules": len(self._module_cache) - sum(self._module_cache.values())
        }


class ModuleFileProcessor:
    """Processor for handling module file operations."""
    
    def __init__(self, logger):
        """Initialize file processor."""
        self.logger = logger
    
    def extract_module_id_from_filename(self, filename: str) -> Optional[str]:
        """
        Extract module ID from filename.
        
        Args:
            filename: CSV filename
            
        Returns:
            Module ID if found, None otherwise
        """
        try:
            # Handle different filename patterns
            # e.g., "organism_M00001.csv" -> "M00001"
            if "_" in filename:
                parts = filename.split("_")
                for part in parts:
                    if part.startswith("M") and part.replace(".csv", "").replace("M", "").isdigit():
                        return part.replace(".csv", "")
            
            # Handle direct module ID filenames
            # e.g., "M00001.csv" -> "M00001"
            base_name = filename.replace(".csv", "")
            if base_name.startswith("M") and base_name[1:].isdigit():
                return base_name
            
            return None
            
        except Exception as e:
            self.logger.warning(f"Failed to extract module ID from filename", 
                              filename=filename, 
                              error=str(e))
            return None
    
    def copy_file_to_target(self, source_path: Path, target_path: Path) -> bool:
        """
        Copy file from source to target with error handling.
        
        Args:
            source_path: Source file path
            target_path: Target file path
            
        Returns:
            True if successful
        """
        try:
            # Ensure target directory exists
            target_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Copy file
            shutil.copy2(source_path, target_path)
            
            self.logger.debug(f"File copied successfully", 
                            source=str(source_path),
                            target=str(target_path))
            
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to copy file", 
                            source=str(source_path),
                            target=str(target_path),
                            error=str(e))
            return False


class KEGGPathwayExtractor(FileSystemProcessor):
    """
    Professional KEGG pathway extractor with metabolic module filtering.
    """
    
    def __init__(self, 
                 config: BioinformaticsConfig,
                 extraction_config: ExtractionConfig):
        """
        Initialize KEGG pathway extractor.
        
        Args:
            config: Base processing configuration
            extraction_config: Extraction-specific configuration
        """
        super().__init__(config, name="KEGGPathwayExtractor")
        
        self.extraction_config = extraction_config
        self.api_client = KEGGModuleAPIClient(
            rate_limit_delay=extraction_config.rate_limit_delay,
            max_retries=extraction_config.max_retries,
            timeout=extraction_config.timeout
        )
        self.file_processor = ModuleFileProcessor(self.logger)
        
        # Statistics tracking
        self.extraction_stats = {
            "files_found": 0,
            "files_processed": 0,
            "metabolic_modules_found": 0,
            "files_copied": 0,
            "api_requests_made": 0,
            "errors_encountered": 0,
        }
    
    def process(self, input_data: Any = None, **kwargs) -> ProcessingResult:
        """
        Extract metabolic modules from KEGG pathway data.
        
        Args:
            input_data: Optional input data (not used in this implementation)
            **kwargs: Additional processing parameters
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            self.logger.info(f"Starting KEGG pathway extraction",
                           source=str(self.extraction_config.source_folder),
                           target=str(self.extraction_config.target_folder))
            
            # Ensure target directory exists
            self.extraction_config.target_folder.mkdir(parents=True, exist_ok=True)
            
            # Find all CSV files in source directory
            csv_files = self._find_csv_files()
            self.extraction_stats["files_found"] = len(csv_files)
            
            if not csv_files:
                result.add_warning("No CSV files found in source directory")
                return result
            
            self.logger.info(f"Found CSV files to process", count=len(csv_files))
            
            # Process files in parallel
            processing_results = self._process_files_parallel(csv_files)
            
            # Aggregate results
            successful_extractions = 0
            metabolic_modules = []
            
            for file_result in processing_results:
                if file_result.success:
                    successful_extractions += 1
                    if file_result.metadata.get("is_metabolic", False):
                        metabolic_modules.append(file_result.metadata.get("module_info"))
                else:
                    result.errors.extend(file_result.errors)
                    result.warnings.extend(file_result.warnings)
            
            # Update statistics
            self.extraction_stats["files_processed"] = successful_extractions
            self.extraction_stats["metabolic_modules_found"] = len(metabolic_modules)
            
            result.metadata.update({
                "total_files_found": len(csv_files),
                "files_processed": successful_extractions,
                "metabolic_modules_found": len(metabolic_modules),
                "metabolic_modules": metabolic_modules,
                "extraction_stats": self.extraction_stats,
                "api_cache_stats": self.api_client.get_cache_stats()
            })
            
            if successful_extractions == 0:
                result.add_error("No files were processed successfully")
            
            self.logger.info(f"KEGG pathway extraction completed",
                           processed=successful_extractions,
                           metabolic=len(metabolic_modules))
            
        except Exception as e:
            result.add_error(f"Extraction failed: {str(e)}")
            self.extraction_stats["errors_encountered"] += 1
        
        return result
    
    def _find_csv_files(self) -> List[Path]:
        """
        Find all CSV files in source directory.
        
        Returns:
            List of CSV file paths
        """
        try:
            csv_files = []
            
            for root, dirs, files in os.walk(self.extraction_config.source_folder):
                for file in files:
                    if (file.endswith(".csv") and 
                        not file.startswith("plants_list")):
                        csv_files.append(Path(root) / file)
            
            self.logger.debug(f"Found CSV files", count=len(csv_files))
            return csv_files
            
        except Exception as e:
            self.logger.error(f"Failed to find CSV files", error=str(e))
            raise
    
    def _process_files_parallel(self, csv_files: List[Path]) -> List[ProcessingResult]:
        """
        Process files in parallel.
        
        Args:
            csv_files: List of CSV files to process
            
        Returns:
            List of processing results
        """
        if self.config.max_workers == 1:
            # Sequential processing
            return [self._process_single_file(file_path) for file_path in csv_files]
        
        # Parallel processing
        results = []
        
        with ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
            future_to_file = {
                executor.submit(self._process_single_file, file_path): file_path
                for file_path in csv_files
            }
            
            for future in as_completed(future_to_file):
                file_path = future_to_file[future]
                try:
                    result = future.result(timeout=self.config.timeout)
                    results.append(result)
                except Exception as e:
                    error_result = ProcessingResult(success=False)
                    error_result.add_error(f"Processing failed for {file_path}: {str(e)}")
                    results.append(error_result)
        
        return results
    
    def _process_single_file(self, file_path: Path) -> ProcessingResult:
        """
        Process single CSV file.
        
        Args:
            file_path: Path to CSV file
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            # Extract module ID from filename
            module_id = self.file_processor.extract_module_id_from_filename(file_path.name)
            
            if not module_id:
                result.add_warning(f"Could not extract module ID from filename: {file_path.name}")
                return result
            
            # Check if module is metabolic
            is_metabolic = False
            if self.extraction_config.verify_metabolic:
                is_metabolic = self.api_client.is_metabolic_module(module_id)
                self.extraction_stats["api_requests_made"] += 1
            else:
                # Assume all modules are metabolic if verification disabled
                is_metabolic = True
            
            module_info = MetabolicModuleInfo(
                module_id=module_id,
                is_metabolic=is_metabolic,
                file_path=file_path
            )
            
            # Copy file if it's metabolic
            if is_metabolic and self.extraction_config.copy_files:
                # Preserve directory structure in target
                relative_path = file_path.relative_to(self.extraction_config.source_folder)
                target_path = self.extraction_config.target_folder / relative_path
                
                if self.file_processor.copy_file_to_target(file_path, target_path):
                    self.extraction_stats["files_copied"] += 1
                    module_info.file_path = target_path
                    
                    self.logger.debug(f"Metabolic module file copied",
                                    module=module_id,
                                    source=str(file_path),
                                    target=str(target_path))
                else:
                    result.add_error(f"Failed to copy metabolic module file: {file_path}")
            
            result.metadata.update({
                "module_id": module_id,
                "is_metabolic": is_metabolic,
                "module_info": module_info,
                "source_path": str(file_path)
            })
            
            if is_metabolic:
                self.logger.debug(f"Metabolic module identified", module=module_id)
            
        except Exception as e:
            result.add_error(f"File processing failed: {str(e)}")
            self.extraction_stats["errors_encountered"] += 1
        
        return result
    
    def get_extraction_statistics(self) -> Dict[str, Any]:
        """Get comprehensive extraction statistics."""
        base_stats = self.get_statistics()
        api_stats = self.api_client.get_cache_stats()
        
        return {
            **base_stats,
            **self.extraction_stats,
            "api_cache_stats": api_stats,
            "extraction_config": {
                "source_folder": str(self.extraction_config.source_folder),
                "target_folder": str(self.extraction_config.target_folder),
                "verify_metabolic": self.extraction_config.verify_metabolic,
                "copy_files": self.extraction_config.copy_files,
            }
        }
    
    def validate_input(self, input_data: Any = None) -> bool:
        """
        Validate extraction configuration.
        
        Args:
            input_data: Input data (not used)
            
        Returns:
            True if configuration is valid
        """
        return (self.extraction_config.source_folder.exists() and
                self.extraction_config.source_folder.is_dir())


def main():
    """Main entry point for KEGG pathway extraction."""
    # Configuration
    source_folder = Path("/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules")
    target_folder = Path("/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic")
    
    # Setup configuration
    config = BioinformaticsConfig(
        input_dir=source_folder,
        output_dir=target_folder,
        max_workers=10,
        log_level="INFO"
    )
    
    extraction_config = ExtractionConfig(
        source_folder=source_folder,
        target_folder=target_folder,
        rate_limit_delay=0.35,
        copy_files=True,
        verify_metabolic=True
    )
    
    # Initialize extractor
    try:
        extractor = KEGGPathwayExtractor(
            config=config,
            extraction_config=extraction_config
        )
        
        # Run extraction
        result = extractor.run()
        
        # Display results
        if result.success:
            print("✅ KEGG pathway extraction completed successfully!")
            print(f"   Files found: {result.metadata.get('total_files_found', 0)}")
            print(f"   Files processed: {result.metadata.get('files_processed', 0)}")
            print(f"   Metabolic modules found: {result.metadata.get('metabolic_modules_found', 0)}")
            print(f"   Processing time: {result.processing_time:.2f} seconds")
        else:
            print("❌ KEGG pathway extraction failed!")
            for error in result.errors:
                print(f"   Error: {error}")
        
        # Display warnings
        for warning in result.warnings:
            print(f"   Warning: {warning}")
        
        # Display statistics
        stats = extractor.get_extraction_statistics()
        print(f"\nExtraction Statistics:")
        print(f"   API requests made: {stats['api_requests_made']}")
        print(f"   Files copied: {stats['files_copied']}")
        print(f"   Cache size: {stats['api_cache_stats']['cache_size']}")
        print(f"   Metabolic modules in cache: {stats['api_cache_stats']['metabolic_modules']}")
        
    except Exception as e:
        print(f"❌ Critical error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())