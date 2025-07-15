"""
Configuration settings for Plant MGC Analysis Pipeline.

This module provides centralized configuration management using pydantic-settings
for environment variables, file-based configuration, and defaults.
"""

import os
from pathlib import Path
from typing import Optional, Dict, Any, List
from functools import lru_cache

from pydantic import BaseModel, Field, validator
from pydantic_settings import BaseSettings, SettingsConfigDict


class DatabaseSettings(BaseModel):
    """Database connection settings."""
    
    # KEGG database settings
    kegg_api_url: str = "https://rest.kegg.jp"
    kegg_request_delay: float = 0.1  # Seconds between requests
    kegg_max_retries: int = 3
    
    # MiBIG database settings  
    mibig_database_path: Optional[Path] = None
    mibig_blast_db: Optional[Path] = None
    mibig_prot_seqs: Optional[Path] = None
    
    # PlantCyc settings
    plantcyc_database_path: Optional[Path] = None
    
    # Ensembl settings
    ensembl_api_url: str = "https://rest.ensembl.org"
    ensembl_request_delay: float = 0.1
    
    # Phytozome settings
    phytozome_api_url: str = "https://phytozome-next.jgi.doe.gov"
    phytozome_token: Optional[str] = None
    
    # NCBI settings
    ncbi_api_url: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    ncbi_email: Optional[str] = None
    ncbi_api_key: Optional[str] = None


class AnalysisSettings(BaseModel):
    """Analysis pipeline settings."""
    
    # Sliding window analysis
    default_window_size: int = 50000  # Base pairs
    default_step_size: int = 10000  # Base pairs
    min_genes_per_window: int = 3
    
    # BLAST settings
    blast_evalue: float = 1e-5
    blast_max_target_seqs: int = 100
    blast_num_threads: int = 4
    
    # Statistical analysis
    significance_threshold: float = 0.05
    multiple_testing_correction: str = "fdr_bh"  # benjamini-hochberg
    
    # Machine learning
    ml_model_path: Optional[Path] = None
    ml_feature_selection: bool = True
    ml_cross_validation_folds: int = 5
    
    # Phylogenetic analysis
    phylo_bootstrap_replicates: int = 1000
    phylo_substitution_model: str = "JTT"


class ComputeSettings(BaseModel):
    """Computational resource settings."""
    
    # Parallel processing
    max_workers: int = Field(default_factory=lambda: os.cpu_count() or 4)
    chunk_size: int = 1000
    
    # Memory management
    max_memory_gb: float = 8.0
    temp_dir: Path = Field(default_factory=lambda: Path("/tmp"))
    
    # SLURM cluster settings
    slurm_partition: str = "general"
    slurm_time_limit: str = "24:00:00"
    slurm_memory_per_cpu: str = "4G"
    use_slurm: bool = False


class LoggingSettings(BaseModel):
    """Logging configuration."""
    
    level: str = "INFO"
    format: str = "{time:YYYY-MM-DD HH:mm:ss} | {level} | {name}:{function}:{line} | {message}"
    log_file: Optional[Path] = None
    max_log_size: str = "10MB"
    log_rotation: str = "7 days"
    
    # Structured logging
    enable_json_logging: bool = False
    log_to_console: bool = True
    log_to_file: bool = True


class DataPathSettings(BaseModel):
    """Data path settings - replaces hardcoded paths from 84 legacy scripts."""
    
    # Base directories (replaces /groups/itay_mayrose/alongonda/...)
    base_data_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda"))
    desktop_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/desktop"))
    datasets_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/datasets"))
    
    # Genome data directories
    full_genomes_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/datasets/full_genomes"))
    ensembl_genomes_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/datasets/full_genomes/ensembl"))
    phytozome_genomes_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/datasets/full_genomes/phytozome"))
    ncbi_genomes_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/datasets/full_genomes/ncbi"))
    
    # Specific organism directories
    arabidopsis_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/desktop/arabidopsis"))
    rice_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/desktop/rice"))
    maize_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/desktop/maize"))
    
    # MGC and reference data
    mgc_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/desktop/MGCs"))
    mibig_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/desktop/MGCs/all_genes_from_mibig"))
    plant_mgc_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/Plant_MGC"))
    
    # Output directories
    results_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/desktop/results"))
    blast_results_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/desktop/blast_results"))
    sliding_window_output_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/desktop/sliding_window_output"))
    
    # Temporary and cache directories
    temp_dir: Path = Field(default_factory=lambda: Path("/tmp/plant_mgc_analysis"))
    cache_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/desktop/.cache"))
    
    # Configuration files
    config_dir: Path = Field(default_factory=lambda: Path("/groups/itay_mayrose/alongonda/desktop/config"))
    
    def get_organism_dir(self, organism: str) -> Path:
        """Get organism-specific directory."""
        organism_lower = organism.lower().replace(" ", "_")
        return self.desktop_dir / organism_lower
    
    def get_genome_file_path(self, organism: str, file_type: str = "genome") -> Path:
        """Get genome file path for organism."""
        org_dir = self.get_organism_dir(organism)
        
        if file_type == "genome":
            return org_dir / f"{organism}.genome.fasta"
        elif file_type == "protein":
            return org_dir / f"{organism}.protein.fasta"
        elif file_type == "annotation":
            return org_dir / f"{organism}.gff3"
        else:
            raise ValueError(f"Unknown file type: {file_type}")
    
    def get_analysis_output_dir(self, analysis_type: str, organism: str = None) -> Path:
        """Get analysis-specific output directory."""
        base_dir = self.results_dir / analysis_type
        if organism:
            base_dir = base_dir / organism
        return base_dir


class Settings(BaseSettings):
    """Main application settings."""
    
    model_config = SettingsConfigDict(
        env_file=".env",
        env_prefix="PLANT_MGC_",
        case_sensitive=False,
        env_nested_delimiter="__",
        extra="allow",
    )
    
    # Application metadata
    app_name: str = "Plant MGC Analysis Pipeline"
    app_version: str = "0.1.0"
    debug: bool = False
    
    # Data directories
    data_dir: Path = Field(default_factory=lambda: Path.cwd() / "data")
    output_dir: Path = Field(default_factory=lambda: Path.cwd() / "output")
    cache_dir: Path = Field(default_factory=lambda: Path.cwd() / ".cache")
    
    # Data paths (replaces hardcoded paths)
    paths: DataPathSettings = Field(default_factory=DataPathSettings)
    
    # Configuration sections
    database: DatabaseSettings = Field(default_factory=DatabaseSettings)
    analysis: AnalysisSettings = Field(default_factory=AnalysisSettings)
    compute: ComputeSettings = Field(default_factory=ComputeSettings)
    logging: LoggingSettings = Field(default_factory=LoggingSettings)
    
    # Custom configuration
    custom_config: Dict[str, Any] = Field(default_factory=dict)
    
    # Environment detection
    is_cluster: bool = Field(default_factory=lambda: Path("/groups").exists())
    is_local: bool = Field(default_factory=lambda: not Path("/groups").exists())
    
    @validator("data_dir", "output_dir", "cache_dir", pre=True)
    def create_directories(cls, v):
        """Create directories if they don't exist."""
        path = Path(v)
        path.mkdir(parents=True, exist_ok=True)
        return path
    
    @validator("paths", pre=True)
    def setup_data_paths(cls, v):
        """Setup data paths based on environment."""
        if isinstance(v, dict):
            return DataPathSettings(**v)
        elif isinstance(v, DataPathSettings):
            return v
        else:
            return DataPathSettings()
    
    @validator("database", pre=True)
    def validate_database_paths(cls, v):
        """Validate database paths exist if provided."""
        if isinstance(v, dict):
            for key, value in v.items():
                if key.endswith("_path") and value is not None:
                    path = Path(value)
                    if not path.exists():
                        raise ValueError(f"Database path does not exist: {path}")
        return v
    
    def get_database_config(self) -> Dict[str, Any]:
        """Get database configuration as dictionary."""
        return self.database.model_dump()
    
    def get_analysis_config(self) -> Dict[str, Any]:
        """Get analysis configuration as dictionary."""
        return self.analysis.model_dump()
    
    def get_compute_config(self) -> Dict[str, Any]:
        """Get compute configuration as dictionary."""
        return self.compute.model_dump()
    
    def get_logging_config(self) -> Dict[str, Any]:
        """Get logging configuration as dictionary."""
        return self.logging.model_dump()
    
    def get_paths_config(self) -> Dict[str, Any]:
        """Get paths configuration as dictionary."""
        return self.paths.model_dump()
    
    def get_organism_config(self, organism: str) -> Dict[str, Any]:
        """Get organism-specific configuration."""
        return {
            "organism": organism,
            "organism_dir": str(self.paths.get_organism_dir(organism)),
            "genome_file": str(self.paths.get_genome_file_path(organism, "genome")),
            "protein_file": str(self.paths.get_genome_file_path(organism, "protein")),
            "annotation_file": str(self.paths.get_genome_file_path(organism, "annotation")),
            "output_dir": str(self.paths.get_analysis_output_dir("general", organism)),
        }
    
    def migrate_legacy_paths(self) -> Dict[str, str]:
        """Get mapping of legacy hardcoded paths to new configuration."""
        return {
            "/groups/itay_mayrose/alongonda/datasets/": str(self.paths.datasets_dir),
            "/groups/itay_mayrose/alongonda/desktop/": str(self.paths.desktop_dir),
            "/groups/itay_mayrose/alongonda/Plant_MGC/": str(self.paths.plant_mgc_dir),
            "/groups/itay_mayrose/alongonda/desktop/MGCs/": str(self.paths.mgc_dir),
            "/groups/itay_mayrose/alongonda/desktop/arabidopsis/": str(self.paths.arabidopsis_dir),
            "/groups/itay_mayrose/alongonda/desktop/results/": str(self.paths.results_dir),
        }
    
    def update_config(self, **updates: Any) -> None:
        """Update configuration with new values."""
        for key, value in updates.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                self.custom_config[key] = value
    
    def save_config(self, path: Path) -> None:
        """Save current configuration to file."""
        import json
        with open(path, 'w') as f:
            json.dump(self.model_dump(), f, indent=2, default=str)
    
    @classmethod
    def load_config(cls, path: Path) -> "Settings":
        """Load configuration from file."""
        import json
        with open(path, 'r') as f:
            config_data = json.load(f)
        return cls(**config_data)


@lru_cache()
def get_settings() -> Settings:
    """Get cached settings instance."""
    return Settings()


def get_config_template() -> Dict[str, Any]:
    """Get a template configuration dictionary."""
    return Settings().model_dump()


def create_environment_config(environment: str = "production") -> Settings:
    """Create environment-specific configuration."""
    if environment == "development":
        return Settings(
            debug=True,
            logging=LoggingSettings(level="DEBUG"),
            compute=ComputeSettings(max_workers=2),
        )
    elif environment == "cluster":
        return Settings(
            compute=ComputeSettings(
                use_slurm=True,
                max_workers=16,
                slurm_partition="general",
                slurm_time_limit="48:00:00",
                slurm_memory_per_cpu="8G",
            ),
        )
    else:
        return Settings()


def migrate_legacy_script_paths(script_content: str, settings: Settings) -> str:
    """Migrate hardcoded paths in legacy script content."""
    path_mapping = settings.migrate_legacy_paths()
    
    for old_path, new_path in path_mapping.items():
        script_content = script_content.replace(old_path, new_path)
    
    return script_content