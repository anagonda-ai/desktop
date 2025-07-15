#!/usr/bin/env python3
"""
State-of-the-art KEGG metabolic pathway extractor for local processing.

This module processes KEGG annotation files to extract metabolic modules with
professional error handling, caching, and comprehensive validation.

The script scans a directory tree for CSV files containing KEGG module data,
checks if each module is metabolic, and copies metabolic modules to a separate
directory structure.

Usage:
    python kegg_pathway_extractor_local.py [source_folder] [output_folder]

Example:
    python kegg_pathway_extractor_local.py /data/kegg_modules /data/metabolic_modules
"""

import sys
import os
from pathlib import Path
from typing import Dict, List, Optional, Set, Any
from dataclasses import dataclass
import shutil
import logging

# Import our new framework
sys.path.insert(0, str(Path(__file__).parent.parent))
from core import BaseProcessor, Config, get_config
from core.exceptions import ProcessingError, ValidationError, FileOperationError
from core.types import ProcessingConfig, AnalysisResult, ProcessingStatus
from utils import KEGGClient, FileManager, BatchProcessor, ProgressTracker, CacheManager


@dataclass
class MetabolicExtractionConfig(ProcessingConfig):
    """Configuration for metabolic module extraction."""
    max_workers: int = 10
    cache_metabolic_results: bool = True
    copy_files: bool = True
    validate_csv_format: bool = True
    exclude_patterns: List[str] = None
    
    def __post_init__(self):
        if self.exclude_patterns is None:
            self.exclude_patterns = ["plants_list"]


@dataclass
class ModuleInfo:
    """Information about a KEGG module."""
    module_id: str
    file_path: Path
    is_metabolic: Optional[bool] = None
    class_info: Optional[str] = None
    
    @property
    def target_path(self) -> Path:
        """Get target path for metabolic module."""
        return Path(str(self.file_path).replace(
            "KEGG_annotations_modules", 
            "KEGG_annotations_modules_metabolic"
        ))


class MetabolicModuleChecker:
    """
    Professional metabolic module classification system.
    
    Checks KEGG modules to determine if they are metabolic using
    caching and batch processing for efficiency.
    """
    
    def __init__(self, config: Optional[MetabolicExtractionConfig] = None):
        self.config = config or MetabolicExtractionConfig()
        self.kegg_client = KEGGClient()
        self.cache = CacheManager(
            cache_dir=Path("./cache/metabolic_modules"),
            ttl_hours=24 * 7  # Cache for a week
        )
        self.logger = logging.getLogger(__name__)
        
        # In-memory cache for current session
        self.metabolic_results: Dict[str, bool] = {}
    
    def is_metabolic_module(self, module_id: str) -> bool:
        """
        Check if a module is metabolic with caching.
        
        Args:
            module_id: KEGG module ID (e.g., 'M00001')
            
        Returns:
            True if module is metabolic, False otherwise
        """
        # Check session cache first
        if module_id in self.metabolic_results:
            self.logger.debug(f"Session cache hit for module {module_id}")
            return self.metabolic_results[module_id]
        
        # Check persistent cache
        cache_key = f"metabolic:{module_id}"
        cached_result = self.cache.get(cache_key)
        if cached_result is not None:
            self.logger.debug(f"Persistent cache hit for module {module_id}")
            self.metabolic_results[module_id] = cached_result
            return cached_result
        
        # Fetch from KEGG API
        try:
            self.logger.debug(f"Fetching module classification for {module_id}")
            response = self.kegg_client.get(f"get/{module_id}")
            
            if response.status_code != 200:
                self.logger.warning(f"Failed to fetch module {module_id}: status {response.status_code}")
                self.metabolic_results[module_id] = False
                return False
            
            # Parse CLASS field to determine if metabolic
            is_metabolic = self._parse_metabolic_class(response.data)
            
            # Cache results
            self.metabolic_results[module_id] = is_metabolic
            self.cache.set(cache_key, is_metabolic)
            
            self.logger.debug(f"Module {module_id} is {'metabolic' if is_metabolic else 'non-metabolic'}")
            return is_metabolic
            
        except Exception as e:
            self.logger.error(f"Error checking module {module_id}: {e}")
            self.metabolic_results[module_id] = False
            return False
    
    def _parse_metabolic_class(self, kegg_data: str) -> bool:
        """
        Parse KEGG module data to determine if metabolic.
        
        Args:
            kegg_data: Raw KEGG module data
            
        Returns:
            True if module contains metabolic keywords
        """
        try:
            for line in kegg_data.splitlines():
                if line.startswith("CLASS"):
                    class_info = line[12:].strip().lower()
                    # Check for metabolic keywords
                    metabolic_keywords = ["metabol", "biosynthesis", "degradation", "pathway"]
                    if any(keyword in class_info for keyword in metabolic_keywords):
                        return True
            return False
            
        except Exception as e:
            self.logger.error(f"Error parsing KEGG class data: {e}")
            return False
    
    def get_module_class_info(self, module_id: str) -> Optional[str]:
        """Get the class information for a module."""
        try:
            response = self.kegg_client.get(f"get/{module_id}")
            if response.status_code == 200:
                for line in response.data.splitlines():
                    if line.startswith("CLASS"):
                        return line[12:].strip()
            return None
        except Exception:
            return None


class KEGGMetabolicExtractor(BaseProcessor):
    """
    State-of-the-art KEGG metabolic module extractor.
    
    Processes KEGG annotation files to identify and extract metabolic modules
    with comprehensive error handling, progress tracking, and validation.
    """
    
    def __init__(self, config: Optional[MetabolicExtractionConfig] = None):
        super().__init__(config)
        self.config = config or MetabolicExtractionConfig()
        
        # Initialize components
        self.file_manager = FileManager()
        self.metabolic_checker = MetabolicModuleChecker(self.config)
        self.batch_processor = BatchProcessor(
            max_workers=self.config.max_workers,
            executor_type="thread"  # Better for I/O operations
        )
        
        # Statistics
        self.processed_files = 0
        self.metabolic_modules_found = 0
        self.files_copied = 0
    
    def validate_input(self, data: Dict[str, Any]) -> None:
        """
        Validate input parameters.
        
        Args:
            data: Input parameters dictionary
            
        Raises:
            ValidationError: If validation fails
        """
        # Validate source folder
        source_folder = Path(data.get("source_folder", ""))
        if not source_folder.exists():
            raise ValidationError(
                f"Source folder does not exist: {source_folder}",
                field_name="source_folder",
                invalid_value=str(source_folder)
            )
        
        if not source_folder.is_dir():
            raise ValidationError(
                f"Source path is not a directory: {source_folder}",
                field_name="source_folder"
            )
        
        # Output folder will be created if it doesn't exist
        output_folder = Path(data.get("output_folder", ""))
        if output_folder.exists() and not output_folder.is_dir():
            raise ValidationError(
                f"Output path exists but is not a directory: {output_folder}",
                field_name="output_folder"
            )
    
    def discover_csv_files(self, source_folder: Path) -> List[Path]:
        """
        Discover all CSV files in the source directory tree.
        
        Args:
            source_folder: Root folder to scan
            
        Returns:
            List of CSV file paths
        """
        self.logger.info(f"Scanning for CSV files in: {source_folder}")
        
        csv_files = []
        exclude_patterns = self.config.exclude_patterns
        
        try:
            # Use Path.rglob for recursive file discovery
            for csv_file in source_folder.rglob("*.csv"):
                # Check exclusion patterns
                if any(pattern in csv_file.name for pattern in exclude_patterns):
                    self.logger.debug(f"Excluding file: {csv_file.name}")
                    continue
                
                csv_files.append(csv_file)
            
            self.logger.info(f"Found {len(csv_files)} CSV files to process")
            return csv_files
            
        except Exception as e:
            raise ProcessingError(
                f"Failed to discover CSV files: {e}",
                processor_name=self.__class__.__name__,
                stage="file_discovery",
                original_error=e
            )
    
    def extract_module_id(self, csv_file: Path) -> Optional[str]:
        """
        Extract module ID from CSV filename.
        
        Args:
            csv_file: Path to CSV file
            
        Returns:
            Module ID if found, None otherwise
        """
        try:
            # Expected format: organism_module.csv (e.g., "hsa_M00001.csv")
            filename = csv_file.stem
            parts = filename.split("_")
            
            if len(parts) >= 2:
                module_part = parts[1]
                # Module IDs typically start with 'M'
                if module_part.startswith('M') and len(module_part) >= 6:
                    return module_part
            
            self.logger.warning(f"Could not extract module ID from filename: {filename}")
            return None
            
        except Exception as e:
            self.logger.error(f"Error extracting module ID from {csv_file}: {e}")
            return None
    
    def validate_csv_file(self, csv_file: Path) -> bool:
        """
        Validate CSV file format and content.
        
        Args:
            csv_file: Path to CSV file
            
        Returns:
            True if file is valid
        """
        if not self.config.validate_csv_format:
            return True
        
        try:
            # Basic validation - check if file is readable and not empty
            if csv_file.stat().st_size == 0:
                self.logger.warning(f"Empty CSV file: {csv_file}")
                return False
            
            # Try to read first line to check format
            with open(csv_file, 'r') as f:
                first_line = f.readline().strip()
                if not first_line:
                    self.logger.warning(f"CSV file has no content: {csv_file}")
                    return False
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error validating CSV file {csv_file}: {e}")
            return False
    
    def process_single_file(self, csv_file: Path, output_folder: Path) -> ModuleInfo:
        """
        Process a single CSV file to check if it contains a metabolic module.
        
        Args:
            csv_file: Path to CSV file
            output_folder: Output folder for metabolic modules
            
        Returns:
            ModuleInfo object with processing results
        """
        # Extract module ID
        module_id = self.extract_module_id(csv_file)
        module_info = ModuleInfo(
            module_id=module_id or "unknown",
            file_path=csv_file
        )
        
        if not module_id:
            self.logger.warning(f"Skipping file with no module ID: {csv_file}")
            return module_info
        
        # Validate CSV file
        if not self.validate_csv_file(csv_file):
            self.logger.warning(f"Skipping invalid CSV file: {csv_file}")
            return module_info
        
        try:
            # Check if module is metabolic
            is_metabolic = self.metabolic_checker.is_metabolic_module(module_id)
            module_info.is_metabolic = is_metabolic
            
            if is_metabolic:
                # Get class information for logging
                module_info.class_info = self.metabolic_checker.get_module_class_info(module_id)
                
                if self.config.copy_files:
                    # Create target directory structure
                    target_path = self._get_target_path(csv_file, output_folder)
                    target_path.parent.mkdir(parents=True, exist_ok=True)
                    
                    # Copy file
                    shutil.copy2(csv_file, target_path)
                    module_info.target_path = target_path
                    
                    self.logger.info(f"✅ Metabolic module found: {module_id} -> {target_path}")
                else:
                    self.logger.info(f"✅ Metabolic module identified: {module_id}")
            
            return module_info
            
        except Exception as e:
            self.logger.error(f"Error processing file {csv_file}: {e}")
            return module_info
    
    def _get_target_path(self, source_file: Path, output_folder: Path) -> Path:
        """
        Calculate target path for copied file.
        
        Args:
            source_file: Original file path
            output_folder: Output folder
            
        Returns:
            Target file path
        """
        # Preserve directory structure relative to source
        # Convert absolute path to relative structure
        source_str = str(source_file)
        
        # Replace the base directory name to create metabolic version
        if "KEGG_annotations_modules" in source_str:
            target_str = source_str.replace(
                "KEGG_annotations_modules",
                "KEGG_annotations_modules_metabolic"
            )
            return Path(target_str)
        else:
            # Fallback: use output folder + relative path
            return output_folder / source_file.name
    
    def process_file_batch(self, file_batch: List[Path], output_folder: Path) -> List[ModuleInfo]:
        """
        Process a batch of CSV files.
        
        Args:
            file_batch: List of CSV files to process
            output_folder: Output folder for metabolic modules
            
        Returns:
            List of ModuleInfo objects
        """
        results = []
        for csv_file in file_batch:
            try:
                result = self.process_single_file(csv_file, output_folder)
                results.append(result)
            except Exception as e:
                self.logger.error(f"Error in batch processing file {csv_file}: {e}")
                # Add error result
                error_result = ModuleInfo(
                    module_id="error",
                    file_path=csv_file,
                    is_metabolic=False
                )
                results.append(error_result)
        
        return results
    
    def process(self, data: Dict[str, Any], **kwargs) -> Dict[str, Any]:
        """
        Process KEGG annotation files to extract metabolic modules.
        
        Args:
            data: Input parameters
            **kwargs: Additional processing parameters
            
        Returns:
            Dictionary with processing results
        """
        source_folder = Path(data["source_folder"])
        output_folder = Path(data["output_folder"])
        
        # Create output folder
        self.file_manager.ensure_directory(output_folder)
        
        # Discover CSV files
        csv_files = self.discover_csv_files(source_folder)
        
        if not csv_files:
            self.logger.warning("No CSV files found to process")
            return {
                "status": "completed",
                "files_processed": 0,
                "metabolic_modules_found": 0,
                "files_copied": 0
            }
        
        # Process files in batches
        def batch_function(file_batch):
            return self.process_file_batch(file_batch, output_folder)
        
        self.logger.info(f"Processing {len(csv_files)} files with {self.config.max_workers} workers")
        
        batch_results = self.batch_processor.process_parallel(
            csv_files,
            batch_function,
            description="Processing KEGG module files"
        )
        
        # Collect statistics
        total_processed = 0
        metabolic_found = 0
        files_copied = 0
        metabolic_modules = []
        
        for batch_result in batch_results:
            if batch_result.success:
                module_infos = batch_result.result
                for module_info in module_infos:
                    total_processed += 1
                    if module_info.is_metabolic:
                        metabolic_found += 1
                        metabolic_modules.append(module_info)
                        if hasattr(module_info, 'target_path') and module_info.target_path:
                            files_copied += 1
        
        # Log summary
        self.logger.info(
            f"✅ Processing complete: {total_processed} files processed, "
            f"{metabolic_found} metabolic modules found, {files_copied} files copied"
        )
        
        # Log some example metabolic modules
        if metabolic_modules:
            self.logger.info("Sample metabolic modules found:")
            for module_info in metabolic_modules[:5]:  # Show first 5
                class_info = module_info.class_info or "Unknown class"
                self.logger.info(f"  {module_info.module_id}: {class_info}")
        
        return {
            "status": "completed",
            "source_folder": str(source_folder),
            "output_folder": str(output_folder),
            "files_processed": total_processed,
            "metabolic_modules_found": metabolic_found,
            "files_copied": files_copied,
            "metabolic_modules": [
                {
                    "module_id": m.module_id,
                    "file_path": str(m.file_path),
                    "class_info": m.class_info
                }
                for m in metabolic_modules[:10]  # Include first 10 in results
            ]
        }


def main():
    """
    Main entry point for KEGG metabolic extractor.
    
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
        source_folder = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules"
        output_folder = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic"
    elif len(sys.argv) == 3:
        source_folder = sys.argv[1]
        output_folder = sys.argv[2]
    else:
        print("Usage: python kegg_pathway_extractor_local.py [source_folder] [output_folder]")
        print("\nArguments:")
        print("  source_folder  : Directory containing KEGG annotation CSV files")
        print("  output_folder  : Directory for metabolic module output")
        print("\nExample:")
        print("  python kegg_pathway_extractor_local.py /data/kegg_modules /data/metabolic_modules")
        print("\nNote: If no arguments provided, uses default paths for backward compatibility")
        sys.exit(1)
    
    # Create processor with configuration
    config = MetabolicExtractionConfig(
        max_workers=10,
        cache_metabolic_results=True,
        copy_files=True,
        validate_csv_format=True
    )
    
    processor = KEGGMetabolicExtractor(config)
    
    try:
        # Prepare input data
        input_data = {
            "source_folder": source_folder,
            "output_folder": output_folder
        }
        
        # Run processing
        logger.info("Starting KEGG metabolic module extraction")
        result = processor.run(input_data)
        
        if result.status == ProcessingStatus.COMPLETED:
            logger.info("Extraction completed successfully")
            if result.results:
                stats = result.results
                print(f"✅ All metabolic module-based KEGG data collected.")
                print(f"   Files processed: {stats['files_processed']}")
                print(f"   Metabolic modules found: {stats['metabolic_modules_found']}")
                print(f"   Files copied: {stats['files_copied']}")
        else:
            logger.error(f"Extraction failed: {result.error_message}")
            sys.exit(1)
    
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()