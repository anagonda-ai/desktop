"""
Industry-Level Archive Extractor.

This module provides object-oriented archive extraction with comprehensive
error handling, validation, parallel processing, and multiple format support.
"""

import os
import gzip
import shutil
import zipfile
import tarfile
from typing import Dict, List, Optional, Any, Set, Tuple
from dataclasses import dataclass, field
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from enum import Enum

from ..base_classes import (
    BaseProcessor,
    BioinformaticsConfig,
    ProcessingResult,
    FileSystemProcessor,
    ParallelProcessor
)


class ArchiveFormat(Enum):
    """Supported archive formats."""
    GZIP = "gzip"
    ZIP = "zip"
    TAR = "tar"
    TAR_GZ = "tar.gz"
    TAR_BZ2 = "tar.bz2"
    TAR_XZ = "tar.xz"


@dataclass
class ExtractionConfig:
    """Configuration for archive extraction operations."""
    
    supported_formats: List[ArchiveFormat] = field(default_factory=lambda: [
        ArchiveFormat.GZIP, ArchiveFormat.ZIP, ArchiveFormat.TAR,
        ArchiveFormat.TAR_GZ, ArchiveFormat.TAR_BZ2, ArchiveFormat.TAR_XZ
    ])
    delete_after_extraction: bool = False
    preserve_directory_structure: bool = True
    overwrite_existing: bool = False
    verify_extraction: bool = True
    create_extraction_log: bool = True
    max_extraction_size_gb: float = 10.0
    buffer_size: int = 1024 * 1024  # 1MB


@dataclass
class ExtractionOperation:
    """Information about an extraction operation."""
    
    archive_path: Path
    extraction_path: Path
    archive_format: ArchiveFormat
    status: str = "pending"  # pending, success, failed, skipped
    error_message: Optional[str] = None
    extracted_files: List[Path] = field(default_factory=list)
    archive_size_mb: float = 0.0
    extraction_time: float = 0.0
    
    @property
    def archive_name(self) -> str:
        """Get archive filename."""
        return self.archive_path.name
    
    def __post_init__(self):
        """Calculate archive size."""
        if self.archive_path.exists():
            self.archive_size_mb = self.archive_path.stat().st_size / (1024 * 1024)


class ArchiveFormatDetector:
    """Detects archive formats based on file extensions and content."""
    
    def __init__(self, logger):
        """Initialize format detector."""
        self.logger = logger
        self.format_extensions = {
            ArchiveFormat.GZIP: ['.gz'],
            ArchiveFormat.ZIP: ['.zip'],
            ArchiveFormat.TAR: ['.tar'],
            ArchiveFormat.TAR_GZ: ['.tar.gz', '.tgz'],
            ArchiveFormat.TAR_BZ2: ['.tar.bz2', '.tbz2'],
            ArchiveFormat.TAR_XZ: ['.tar.xz', '.txz']
        }
    
    def detect_format(self, file_path: Path) -> Optional[ArchiveFormat]:
        """Detect archive format from file path."""
        try:
            file_name = file_path.name.lower()
            
            # Check for compound extensions first (tar.gz, tar.bz2, etc.)
            for archive_format, extensions in self.format_extensions.items():
                for ext in sorted(extensions, key=len, reverse=True):
                    if file_name.endswith(ext):
                        return archive_format
            
            return None
            
        except Exception as e:
            self.logger.warning(f"Failed to detect archive format", 
                              file=str(file_path), error=str(e))
            return None


class ArchiveExtractor(FileSystemProcessor):
    """
    Professional archive extractor with comprehensive functionality.
    """
    
    def __init__(self, 
                 config: BioinformaticsConfig,
                 extraction_config: ExtractionConfig):
        """
        Initialize archive extractor.
        
        Args:
            config: Base processing configuration
            extraction_config: Extraction-specific configuration
        """
        super().__init__(config, name="ArchiveExtractor")
        
        self.extraction_config = extraction_config
        self.format_detector = ArchiveFormatDetector(self.logger)
        
        # Statistics
        self.extraction_stats = {
            "archives_found": 0,
            "archives_extracted": 0,
            "archives_failed": 0,
            "archives_skipped": 0,
            "total_files_extracted": 0,
            "total_size_extracted_mb": 0.0,
            "archives_deleted": 0,
        }
    
    def process(self, input_data: Path, **kwargs) -> ProcessingResult:
        """
        Process archive extraction operations.
        
        Args:
            input_data: Path to directory containing archives or single archive
            **kwargs: Additional processing parameters including:
                - extraction_dir: Target extraction directory
                - delete_archives: Whether to delete archives after extraction
            
        Returns:
            ProcessingResult object
        """
        result = ProcessingResult(success=True)
        
        try:
            input_path = Path(input_data) if isinstance(input_data, str) else input_data
            extraction_dir = Path(kwargs.get('extraction_dir', input_path))
            delete_archives = kwargs.get('delete_archives', self.extraction_config.delete_after_extraction)
            
            self.logger.info(f"Starting archive extraction", 
                           input=str(input_path),
                           extraction_dir=str(extraction_dir),
                           delete_after=delete_archives)
            
            # Determine if input is file or directory
            if input_path.is_file():
                archive_format = self.format_detector.detect_format(input_path)
                if archive_format:
                    extraction_result = self._extract_single_archive(input_path, extraction_dir, delete_archives)
                    result = result.merge(extraction_result)
                else:
                    result.add_error(f"Unsupported archive format: {input_path}")
            
            elif input_path.is_dir():
                directory_result = self._extract_directory(input_path, extraction_dir, delete_archives)
                result = result.merge(directory_result)
            
            else:
                result.add_error(f"Input path does not exist: {input_path}")
            
            # Update metadata
            result.metadata.update({
                "extraction_stats": self.extraction_stats,
                "input_path": str(input_path),
                "extraction_directory": str(extraction_dir),
                "delete_after_extraction": delete_archives
            })
            
            self.logger.info(f"Archive extraction completed",
                           archives_extracted=self.extraction_stats["archives_extracted"],
                           files_extracted=self.extraction_stats["total_files_extracted"])
            
        except Exception as e:
            result.add_error(f"Archive extraction failed: {str(e)}")
        
        return result
    
    def _extract_directory(self, input_dir: Path, extraction_dir: Path, delete_archives: bool) -> ProcessingResult:
        """Extract all archives in directory."""
        result = ProcessingResult(success=True)
        
        try:
            # Find archive files
            extraction_operations = self._find_archives(input_dir, extraction_dir)
            self.extraction_stats["archives_found"] = len(extraction_operations)
            
            if not extraction_operations:
                result.add_warning("No supported archives found in directory")
                return result
            
            self.logger.info(f"Found archives", count=len(extraction_operations))
            
            # Extract archives in parallel or sequential
            if self.config.max_workers == 1:
                # Sequential processing
                for operation in extraction_operations:
                    extraction_result = self._execute_extraction_operation(operation, delete_archives)
                    if not extraction_result.success:
                        result.errors.extend(extraction_result.errors)
                        result.warnings.extend(extraction_result.warnings)
            else:
                # Parallel processing
                with ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
                    future_to_operation = {
                        executor.submit(self._execute_extraction_operation, operation, delete_archives): operation
                        for operation in extraction_operations
                    }
                    
                    for future in as_completed(future_to_operation):
                        operation = future_to_operation[future]
                        try:
                            extraction_result = future.result(timeout=self.config.timeout)
                            if not extraction_result.success:
                                result.errors.extend(extraction_result.errors)
                                result.warnings.extend(extraction_result.warnings)
                        except Exception as e:
                            result.add_error(f"Failed to extract {operation.archive_name}: {str(e)}")
                            self.extraction_stats["archives_failed"] += 1
            
            return result
            
        except Exception as e:
            result.add_error(f"Directory extraction failed: {str(e)}")
            return result
    
    def _find_archives(self, input_dir: Path, extraction_dir: Path) -> List[ExtractionOperation]:
        """Find and analyze archive files in directory."""
        operations = []
        
        try:
            for root, dirs, files in os.walk(input_dir):
                for filename in files:
                    file_path = Path(root) / filename
                    
                    archive_format = self.format_detector.detect_format(file_path)
                    if archive_format and archive_format in self.extraction_config.supported_formats:
                        
                        # Determine extraction path
                        if self.extraction_config.preserve_directory_structure:
                            relative_path = file_path.relative_to(input_dir)
                            target_dir = extraction_dir / relative_path.parent
                        else:
                            target_dir = extraction_dir
                        
                        operation = ExtractionOperation(
                            archive_path=file_path,
                            extraction_path=target_dir,
                            archive_format=archive_format
                        )
                        
                        # Check size limits
                        if operation.archive_size_mb > self.extraction_config.max_extraction_size_gb * 1024:
                            self.logger.warning(f"Archive exceeds size limit", 
                                              archive=operation.archive_name,
                                              size_mb=operation.archive_size_mb,
                                              limit_gb=self.extraction_config.max_extraction_size_gb)
                            continue
                        
                        operations.append(operation)
                        
                        self.logger.debug(f"Found archive", 
                                        file=operation.archive_name,
                                        format=archive_format.value,
                                        size_mb=operation.archive_size_mb)
            
            return operations
            
        except Exception as e:
            self.logger.error(f"Failed to find archives", error=str(e))
            raise
    
    def _extract_single_archive(self, archive_path: Path, extraction_dir: Path, delete_archive: bool) -> ProcessingResult:
        """Extract single archive file."""
        archive_format = self.format_detector.detect_format(archive_path)
        if not archive_format:
            result = ProcessingResult(success=False)
            result.add_error(f"Could not detect format for {archive_path}")
            return result
        
        operation = ExtractionOperation(
            archive_path=archive_path,
            extraction_path=extraction_dir,
            archive_format=archive_format
        )
        
        return self._execute_extraction_operation(operation, delete_archive)
    
    def _execute_extraction_operation(self, operation: ExtractionOperation, delete_archive: bool) -> ProcessingResult:
        """Execute single extraction operation."""
        result = ProcessingResult(success=True)
        
        try:
            import time
            start_time = time.time()
            
            self.logger.debug(f"Extracting archive", 
                            archive=operation.archive_name,
                            format=operation.archive_format.value)
            
            # Create extraction directory
            operation.extraction_path.mkdir(parents=True, exist_ok=True)
            
            # Perform extraction based on format
            if operation.archive_format == ArchiveFormat.GZIP:
                extracted_files = self._extract_gzip(operation)
            elif operation.archive_format == ArchiveFormat.ZIP:
                extracted_files = self._extract_zip(operation)
            elif operation.archive_format in [ArchiveFormat.TAR, ArchiveFormat.TAR_GZ, 
                                            ArchiveFormat.TAR_BZ2, ArchiveFormat.TAR_XZ]:
                extracted_files = self._extract_tar(operation)
            else:
                raise ValueError(f"Unsupported archive format: {operation.archive_format}")
            
            operation.extracted_files = extracted_files
            operation.extraction_time = time.time() - start_time
            operation.status = "success"
            
            # Verify extraction if configured
            if self.extraction_config.verify_extraction:
                if not self._verify_extraction(operation):
                    result.add_error(f"Extraction verification failed for {operation.archive_name}")
                    operation.status = "failed"
                    return result
            
            # Delete archive if requested
            if delete_archive:
                try:
                    operation.archive_path.unlink()
                    self.extraction_stats["archives_deleted"] += 1
                    self.logger.debug(f"Deleted archive", archive=operation.archive_name)
                except Exception as e:
                    result.add_warning(f"Failed to delete archive {operation.archive_name}: {str(e)}")
            
            # Update statistics
            self.extraction_stats["archives_extracted"] += 1
            self.extraction_stats["total_files_extracted"] += len(extracted_files)
            self.extraction_stats["total_size_extracted_mb"] += operation.archive_size_mb
            
            result.metadata.update({
                "archive_file": str(operation.archive_path),
                "extraction_path": str(operation.extraction_path),
                "files_extracted": len(extracted_files),
                "extraction_time": operation.extraction_time
            })
            
            self.logger.info(f"Archive extracted successfully", 
                           archive=operation.archive_name,
                           files=len(extracted_files),
                           time=f"{operation.extraction_time:.2f}s")
            
        except Exception as e:
            operation.status = "failed"
            operation.error_message = str(e)
            result.add_error(f"Extraction failed for {operation.archive_name}: {str(e)}")
            self.extraction_stats["archives_failed"] += 1
        
        return result
    
    def _extract_gzip(self, operation: ExtractionOperation) -> List[Path]:
        """Extract GZIP archive."""
        output_file = operation.extraction_path / operation.archive_path.stem
        
        with gzip.open(operation.archive_path, 'rb') as f_in:
            with open(output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out, self.extraction_config.buffer_size)
        
        return [output_file]
    
    def _extract_zip(self, operation: ExtractionOperation) -> List[Path]:
        """Extract ZIP archive."""
        extracted_files = []
        
        with zipfile.ZipFile(operation.archive_path, 'r') as zip_ref:
            zip_ref.extractall(operation.extraction_path)
            extracted_files = [operation.extraction_path / name for name in zip_ref.namelist()]
        
        return extracted_files
    
    def _extract_tar(self, operation: ExtractionOperation) -> List[Path]:
        """Extract TAR archive (including compressed variants)."""
        extracted_files = []
        
        # Determine compression mode
        mode_map = {
            ArchiveFormat.TAR: 'r',
            ArchiveFormat.TAR_GZ: 'r:gz',
            ArchiveFormat.TAR_BZ2: 'r:bz2',
            ArchiveFormat.TAR_XZ: 'r:xz'
        }
        
        mode = mode_map.get(operation.archive_format, 'r')
        
        with tarfile.open(operation.archive_path, mode) as tar_ref:
            tar_ref.extractall(operation.extraction_path)
            extracted_files = [operation.extraction_path / member.name for member in tar_ref.getmembers() if member.isfile()]
        
        return extracted_files
    
    def _verify_extraction(self, operation: ExtractionOperation) -> bool:
        """Verify that extraction was successful."""
        try:
            # Check that at least one file was extracted
            if not operation.extracted_files:
                return False
            
            # Check that extracted files exist
            for file_path in operation.extracted_files:
                if not file_path.exists():
                    return False
            
            return True
            
        except Exception:
            return False
    
    def get_extraction_statistics(self) -> Dict[str, Any]:
        """Get comprehensive extraction statistics."""
        base_stats = self.get_statistics()
        
        success_rate = 0.0
        if self.extraction_stats["archives_found"] > 0:
            success_rate = self.extraction_stats["archives_extracted"] / self.extraction_stats["archives_found"]
        
        return {
            **base_stats,
            **self.extraction_stats,
            "success_rate": success_rate,
            "average_files_per_archive": (
                self.extraction_stats["total_files_extracted"] / 
                max(1, self.extraction_stats["archives_extracted"])
            )
        }
    
    def validate_input(self, input_data: Path) -> bool:
        """Validate input data."""
        try:
            input_path = Path(input_data) if isinstance(input_data, str) else input_data
            return input_path.exists()
        except Exception:
            return False


def main():
    """Main entry point for archive extraction."""
    # Configuration
    input_dir = Path("/groups/itay_mayrose/alongonda/datasets/full_genomes/phytozome_test/Phytozome/")
    
    # Setup configuration
    config = BioinformaticsConfig(
        input_dir=input_dir,
        output_dir=input_dir,
        max_workers=4,
        log_level="INFO"
    )
    
    extraction_config = ExtractionConfig(
        delete_after_extraction=True,  # This mimics the original behavior
        preserve_directory_structure=True,
        verify_extraction=True,
        max_extraction_size_gb=10.0
    )
    
    # Initialize extractor
    try:
        extractor = ArchiveExtractor(
            config=config,
            extraction_config=extraction_config
        )
        
        # Run extraction
        result = extractor.run(input_dir, delete_archives=True)
        
        # Display results
        if result.success:
            print("✅ Archive extraction completed successfully!")
            stats = extractor.get_extraction_statistics()
            print(f"   Archives found: {stats['archives_found']}")
            print(f"   Archives extracted: {stats['archives_extracted']}")
            print(f"   Total files extracted: {stats['total_files_extracted']}")
            print(f"   Archives deleted: {stats['archives_deleted']}")
            print(f"   Success rate: {stats['success_rate']:.2%}")
            print(f"   Processing time: {result.processing_time:.2f} seconds")
        else:
            print("❌ Archive extraction failed!")
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