
class HmmerProfileMergerProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: hmmer_profile_merger.py."""
    
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
import gzip
import pyhmmer
import gzip
import shutil
import os

# Directory containing your .hmm files
hmm_dir = "/groups/itay_mayrose/alongonda/datasets/hmmer_profile"  # <-- Change this


for filename in os.listdir(hmm_dir):
    if filename.endswith(".hmm"):
        filepath = os.path.join(hmm_dir, filename)
        
        # Check first two bytes to detect gzip
        with open(filepath, "rb") as f:
            magic = f.read(2)
        
        if magic == b'\x1f\x8b':
            logger.info(f"Decompressing {filename}...")
            with gzip.open(filepath, "rb") as f_in:
                with open(filepath + ".decompressed", "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.replace(filepath + ".decompressed", filepath)
        else:
            logger.info(f"Skipping {filename}: already plain text.")

logger.info("✅ Decompression complete!")

output_dir = os.path.join(hmm_dir, "combined_model")
os.makedirs(output_dir, exist_ok=True)
# Output file for the combined HMM model
output_hmm = os.path.join(output_dir, "combined_model.hmm")

def load_hmms(hmm_dir):
    profiles = []
    for file in os.listdir(hmm_dir):
        if file.endswith((".hmm", ".hmm.gz")):
            path = os.path.join(hmm_dir, file)
            # Open normally or gzip
            open_func = gzip.open if file.endswith(".gz") else open
            with open_func(path, "rb") as f:
                try:
                    for hmm in pyhmmer.plan7.HMMFile(f):
                        profiles.append(hmm)
                except Exception as e:
                    logger.info(f"Skipping {file}: {e}")
    return profiles

def write_combined_hmm(output_path, profiles):
    with open(output_path, "wb") as out_f:
        for hmm in profiles:
            hmm.write(out_f)

hmms = load_hmms(hmm_dir)
write_combined_hmm(output_hmm, hmms)

logger.info(f"Merged {len(hmms)} models into {output_hmm}")

# hmmpress the result
import subprocess
subprocess.run(["hmmpress", output_hmm], check=True)