
class IoUtilsProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: io_utils.py."""
    
    def __init__(self, **kwargs):
        """Initialize processor."""
        super().__init__(**kwargs)
        self.settings = get_settings()
    
    def validate_input(self, data):
        """Validate input data."""
        pass  # Implement validation
    
    def process(self, data, **kwargs):
        """Process data."""
        # Original script logic here
        pass


from ..config.settings import get_settings
from ..core.base import BioinformaticsProcessor
from ..core.exceptions import ProcessingError, ValidationError
from ..utils.file_operations import file_manager
from loguru import logger
import os
import pandas as pd
from joblib import Parallel, delayed

def load_mapping_if_exists(mapping_file):
    if mapping_file and os.path.exists(mapping_file):
        mapping_df = pd.read_csv(mapping_file)
        mapping_df["Organism"] = mapping_df["Organism"].str.replace(" ", "_").str.strip()
        return mapping_df
    return None

def extract_organism_from_subject(subject_gene_str):
    try:
        return subject_gene_str.split(";")[0].split("=")[1].split(".")[0]
    except:
        return "unknown"

def load_one_blast(path):
    if not os.path.exists(path):
        logger.info(f"Warning: file not found: {path}")
        return None
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    df["organism"] = df["subject_gene"].map(extract_organism_from_subject)
    return df

def load_selected_blast_results(comp_df, n_threads=4):
    comp_df = comp_df.dropna(subset=["Directory", "Largest Chromosome File"])
    paths = [os.path.join(str(row["Directory"]), str(row["Largest Chromosome File"]))
             for _, row in comp_df.iterrows()]

    results = Parallel(n_jobs=n_threads)(delayed(load_one_blast)(p) for p in paths)
    results = [df for df in results if df is not None]
    if not results:
        raise ValueError("No valid BLAST result files found.")
    return pd.concat(results, ignore_index=True)
