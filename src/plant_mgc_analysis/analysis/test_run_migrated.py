
class TestRunProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: test_run.py."""
    
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
from cladepp_phylo_profiling_helpfuncs.tree_analysis import analyze_tree_clades_dynamic

TREE_PATH = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/species.nwk"
COMPARISON_CSV = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/without_haaap/output_without_haaap/comparison_results.csv"
MAPPING_FILE = "/groups/itay_mayrose/alongonda/datasets/asaph_aharoni/dataset_organism_mapping.csv"
ANCHOR_GENES = ["adcs", "cs", "pABA_transporter"]
OUTPUT_PREFIX = "/groups/itay_mayrose/alongonda/desktop/tree_test/clade_analysis_test"

os.makedirs(os.path.dirname(OUTPUT_PREFIX), exist_ok=True)

analyze_tree_clades_dynamic(
    tree_path=TREE_PATH,
    comparison_csv=COMPARISON_CSV,
    anchor_genes=ANCHOR_GENES,
    output_prefix=OUTPUT_PREFIX,
    mapping_file=MAPPING_FILE
)

output_file = f"{OUTPUT_PREFIX}_summary.csv"
assert os.path.exists(output_file), "Output file not created"
with open(output_file) as f:
    lines = f.readlines()
    assert len(lines) > 1, "Output file is empty"
logger.info("✅ Test passed")
