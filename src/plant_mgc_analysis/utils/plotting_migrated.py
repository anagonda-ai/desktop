
class PlottingProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: plotting.py."""
    
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
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram

def plot_dendrogram(linkage_matrix, labels):
    plt.figure(figsize=(8, 4))
    dendrogram(linkage_matrix, labels=labels, leaf_rotation=90)
    plt.tight_layout()
    plt.savefig("dendrogram.png")
    plt.close()

def plot_heatmap(npp_matrix):
    sns.clustermap(npp_matrix, cmap="vlag", metric="correlation", method="average", figsize=(10, 8))
    plt.savefig("heatmap.png")
    plt.close()
