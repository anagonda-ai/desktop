"""
Configuration management for bioinformatics operations.

This module provides centralized configuration with support for environment variables,
configuration files, and validation.
"""

import os
from typing import Any, Dict, Optional, List
from dataclasses import dataclass, field
from pathlib import Path
import json
import yaml
from .exceptions import ConfigurationError


@dataclass
class DatabaseConfig:
    """Database connection and API configurations."""
    # KEGG API settings
    kegg_api_url: str = "https://rest.kegg.jp"
    kegg_rate_limit: float = 0.2
    
    # NCBI settings
    ncbi_email: Optional[str] = None
    ncbi_api_key: Optional[str] = None
    ncbi_rate_limit: float = 0.35
    
    # UniProt settings
    uniprot_api_url: str = "https://rest.uniprot.org"
    uniprot_rate_limit: float = 0.1
    
    # Local database paths
    mibig_database_path: Optional[Path] = None
    blast_db_path: Optional[Path] = None
    
    def __post_init__(self):
        # Convert string paths to Path objects
        if self.mibig_database_path and isinstance(self.mibig_database_path, str):
            self.mibig_database_path = Path(self.mibig_database_path)
        if self.blast_db_path and isinstance(self.blast_db_path, str):
            self.blast_db_path = Path(self.blast_db_path)


@dataclass
class ComputeConfig:
    """Computational resource configurations."""
    max_workers: int = 4
    batch_size: int = 100
    memory_limit_gb: Optional[int] = None
    timeout_seconds: int = 3600
    
    # SLURM settings
    use_slurm: bool = False
    slurm_partition: str = "general"
    slurm_job_limit: int = 100
    slurm_memory: str = "8G"
    slurm_time: str = "24:00:00"
    
    # Cache settings
    use_cache: bool = True
    cache_dir: Optional[Path] = None
    cache_ttl_hours: int = 24
    
    def __post_init__(self):
        if self.cache_dir and isinstance(self.cache_dir, str):
            self.cache_dir = Path(self.cache_dir)


@dataclass
class AnalysisConfig:
    """Analysis-specific configurations."""
    # BLAST settings
    blast_evalue: float = 1e-5
    blast_max_targets: int = 100
    blast_word_size: int = 6
    
    # Sliding window settings
    default_window_size: int = 50000
    default_step_size: int = 10000
    min_genes_per_window: int = 3
    
    # Statistical settings
    significance_threshold: float = 0.05
    multiple_testing_correction: str = "fdr_bh"
    
    # MGC detection settings
    mgc_min_genes: int = 3
    mgc_max_gap: int = 20000
    mgc_min_score: float = 0.5


@dataclass
class PathConfig:
    """File path configurations."""
    data_dir: Path = Path("./data")
    output_dir: Path = Path("./output")
    cache_dir: Path = Path("./cache")
    log_dir: Path = Path("./logs")
    temp_dir: Path = Path("./temp")
    
    def __post_init__(self):
        # Ensure all paths are Path objects
        for field_name in ['data_dir', 'output_dir', 'cache_dir', 'log_dir', 'temp_dir']:
            value = getattr(self, field_name)
            if isinstance(value, str):
                setattr(self, field_name, Path(value))
    
    def create_directories(self):
        """Create all configured directories if they don't exist."""
        for field_name in ['data_dir', 'output_dir', 'cache_dir', 'log_dir', 'temp_dir']:
            path = getattr(self, field_name)
            path.mkdir(parents=True, exist_ok=True)


@dataclass
class LoggingConfig:
    """Logging configuration."""
    level: str = "INFO"
    format: str = "{time:YYYY-MM-DD HH:mm:ss} | {level} | {name}:{function}:{line} | {message}"
    rotation: str = "1 week"
    retention: str = "1 month"
    backtrace: bool = True
    diagnose: bool = True
    enable_json: bool = False


class Config:
    """
    Central configuration manager for the bioinformatics toolkit.
    
    Supports loading from environment variables, configuration files,
    and provides validation and defaults.
    """
    
    def __init__(
        self,
        config_file: Optional[Path] = None,
        override_dict: Optional[Dict[str, Any]] = None
    ):
        self.database = DatabaseConfig()
        self.compute = ComputeConfig()
        self.analysis = AnalysisConfig()
        self.paths = PathConfig()
        self.logging = LoggingConfig()
        
        # Load configuration in order of precedence
        self._load_defaults()
        if config_file:
            self._load_from_file(config_file)
        self._load_from_environment()
        if override_dict:
            self._load_from_dict(override_dict)
        
        self._validate_config()
        self.paths.create_directories()
    
    def _load_defaults(self):
        """Load default configuration values."""
        # Default values are already set in dataclass definitions
        pass
    
    def _load_from_file(self, config_file: Path):
        """Load configuration from YAML or JSON file."""
        if not config_file.exists():
            raise ConfigurationError(f"Configuration file not found: {config_file}")
        
        try:
            with open(config_file, 'r') as f:
                if config_file.suffix.lower() in ['.yaml', '.yml']:
                    config_dict = yaml.safe_load(f)
                elif config_file.suffix.lower() == '.json':
                    config_dict = json.load(f)
                else:
                    raise ConfigurationError(f"Unsupported config file format: {config_file.suffix}")
            
            self._load_from_dict(config_dict)
            
        except Exception as e:
            raise ConfigurationError(f"Failed to load configuration from {config_file}: {e}")
    
    def _load_from_environment(self):
        """Load configuration from environment variables."""
        env_mappings = {
            # Database settings
            'BIOINF_KEGG_API_URL': ('database', 'kegg_api_url'),
            'BIOINF_NCBI_EMAIL': ('database', 'ncbi_email'),
            'BIOINF_NCBI_API_KEY': ('database', 'ncbi_api_key'),
            'BIOINF_MIBIG_DB_PATH': ('database', 'mibig_database_path'),
            
            # Compute settings
            'BIOINF_MAX_WORKERS': ('compute', 'max_workers'),
            'BIOINF_BATCH_SIZE': ('compute', 'batch_size'),
            'BIOINF_USE_SLURM': ('compute', 'use_slurm'),
            'BIOINF_SLURM_PARTITION': ('compute', 'slurm_partition'),
            
            # Path settings
            'BIOINF_DATA_DIR': ('paths', 'data_dir'),
            'BIOINF_OUTPUT_DIR': ('paths', 'output_dir'),
            'BIOINF_CACHE_DIR': ('paths', 'cache_dir'),
            
            # Logging settings
            'BIOINF_LOG_LEVEL': ('logging', 'level'),
            'BIOINF_LOG_FORMAT': ('logging', 'format'),
        }
        
        for env_var, (section, field) in env_mappings.items():
            if env_var in os.environ:
                value = os.environ[env_var]
                section_obj = getattr(self, section)
                
                # Type conversion based on field type
                if hasattr(section_obj, field):
                    current_value = getattr(section_obj, field)
                    if isinstance(current_value, bool):
                        value = value.lower() in ('true', '1', 'yes', 'on')
                    elif isinstance(current_value, int):
                        value = int(value)
                    elif isinstance(current_value, float):
                        value = float(value)
                    elif isinstance(current_value, Path):
                        value = Path(value)
                    
                    setattr(section_obj, field, value)
    
    def _load_from_dict(self, config_dict: Dict[str, Any]):
        """Load configuration from dictionary."""
        for section_name, section_config in config_dict.items():
            if hasattr(self, section_name) and isinstance(section_config, dict):
                section_obj = getattr(self, section_name)
                for field_name, value in section_config.items():
                    if hasattr(section_obj, field_name):
                        setattr(section_obj, field_name, value)
    
    def _validate_config(self):
        """Validate configuration values."""
        # Validate required NCBI settings
        if self.database.ncbi_email is None:
            raise ConfigurationError(
                "NCBI email is required. Set BIOINF_NCBI_EMAIL environment variable."
            )
        
        # Validate numeric ranges
        if self.compute.max_workers < 1:
            raise ConfigurationError("max_workers must be at least 1")
        
        if self.compute.batch_size < 1:
            raise ConfigurationError("batch_size must be at least 1")
        
        if not 0 < self.analysis.significance_threshold < 1:
            raise ConfigurationError("significance_threshold must be between 0 and 1")
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            'database': self.database.__dict__,
            'compute': self.compute.__dict__,
            'analysis': self.analysis.__dict__,
            'paths': {k: str(v) for k, v in self.paths.__dict__.items()},
            'logging': self.logging.__dict__
        }
    
    def save_to_file(self, file_path: Path):
        """Save current configuration to file."""
        config_dict = self.to_dict()
        
        with open(file_path, 'w') as f:
            if file_path.suffix.lower() in ['.yaml', '.yml']:
                yaml.dump(config_dict, f, default_flow_style=False)
            elif file_path.suffix.lower() == '.json':
                json.dump(config_dict, f, indent=2)
            else:
                raise ConfigurationError(f"Unsupported file format: {file_path.suffix}")
    
    def load_pipeline_config(self, file_path: Path) -> Dict[str, Any]:
        """Load pipeline configuration from file."""
        if not file_path.exists():
            raise ConfigurationError(f"Pipeline configuration file not found: {file_path}")
        
        try:
            with open(file_path, 'r') as f:
                if file_path.suffix.lower() in ['.yaml', '.yml']:
                    pipeline_config = yaml.safe_load(f)
                elif file_path.suffix.lower() == '.json':
                    pipeline_config = json.load(f)
                else:
                    raise ConfigurationError(f"Unsupported pipeline config format: {file_path.suffix}")
            
            return pipeline_config
            
        except Exception as e:
            raise ConfigurationError(f"Failed to load pipeline configuration: {e}")
    
    def get_workflow_config(self, workflow_name: str) -> Dict[str, Any]:
        """Get workflow-specific configuration."""
        base_config = self.to_dict()
        
        # Add workflow-specific defaults
        workflow_defaults = {
            'kegg': {
                'mode': 'full',
                'organism_batch_size': 10,
                'extract_workers': 10,
                'use_cache': True,
                'skip_existing': True
            },
            'blast': {
                'blast_type': 'blastp',
                'evalue': 1e-3,
                'num_threads': 8,
                'max_target_seqs': 500,
                'min_identity': 30.0,
                'min_coverage': 50.0
            },
            'mgc': {
                'detection_method': 'antiSMASH',
                'min_cluster_size': 3,
                'max_gap_size': 20000
            },
            'phylo': {
                'alignment_method': 'muscle',
                'tree_method': 'fasttree',
                'bootstrap_replicates': 100
            }
        }
        
        if workflow_name.lower() in workflow_defaults:
            base_config['workflow_specific'] = workflow_defaults[workflow_name.lower()]
        
        return base_config


# Global configuration instance
_global_config: Optional[Config] = None


def get_config() -> Config:
    """Get the global configuration instance."""
    global _global_config
    if _global_config is None:
        _global_config = Config()
    return _global_config


def set_config(config: Config):
    """Set the global configuration instance."""
    global _global_config
    _global_config = config


def load_config(config_file: Optional[Path] = None) -> Config:
    """Load and set the global configuration."""
    config = Config(config_file)
    set_config(config)
    return config