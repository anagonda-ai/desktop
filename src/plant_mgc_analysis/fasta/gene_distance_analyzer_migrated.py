
class GeneDistanceAnalyzerProcessor(BioinformaticsProcessor):
    """Migrated from legacy script: gene_distance_analyzer.py."""
    
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
import statistics

def parse_fasta(fasta_file):
    with open(fasta_file, 'r') as file:
        lines = file.readlines()

    mgcs = {}
    current_mgc = None

    for line in lines:
        if line.startswith(">"):
            header = line.strip().split("|")
            mgc_id = header[0][1:]
            chromosome = header[1]
            location = header[2]
            start, end = map(int, location.split("-"))
            strand = header[3]
            gene_id = header[4]
            gene_function = header[5]

            if mgc_id not in mgcs:
                mgcs[mgc_id] = []

            mgcs[mgc_id].append((chromosome, start, end, strand, gene_id, gene_function))
        else:
            continue

    return mgcs

def calculate_max_distance(mgcs):
    max_distances = {}
    max_distance_total = 0
    max_distance_cluster = ""
    distances = []

    for mgc_id, genes in mgcs.items():
        genes.sort(key=lambda x: x[1])  # Sort by start location
        max_distance = 0

        for i in range(len(genes) - 2):
            start1, end1 = genes[i][1], genes[i][2]
            start3, end3 = genes[i+2][1], genes[i+2][2]

            distance = max(start3, end3) - min(start1, end1)
            distances.append(distance)
            
            if distance > max_distance:
                max_distance = distance
            if distance > max_distance_total:
                max_distance_total = distance
                max_distance_cluster = mgc_id

        if not max_distance == 0:
            max_distances[mgc_id] = max_distance
        
    avg_distance = statistics.mean(distances)
    median_distance = statistics.median(distances)
    

    return max_distances, max_distance_total, max_distance_cluster, avg_distance, median_distance

# Path to your FASTA file
fasta_file = input("Enter fasta file path: \n")

mgcs = parse_fasta(fasta_file)
max_distances, max_distance_total, max_distance_cluster, avg_distance, median_distance = calculate_max_distance(mgcs)

for mgc_id, distance in max_distances.items():
    logger.info(f'Maximal distance between 3 neighboring genes in MGC {mgc_id}: {distance}bp')

logger.info(f'Maximal distance between 3 neighboring genes is {max_distance_total}bp in cluster {max_distance_cluster}')
logger.info(f'Average distance between 3 neighboring genes in MGC is {avg_distance}bp')
logger.info(f'Median distance between 3 neighboring genes in MGC is {median_distance}bp')