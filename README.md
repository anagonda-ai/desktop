# Plant MGC Analysis Pipeline

A comprehensive bioinformatics toolkit for identifying and analyzing biosynthetic gene clusters in plant genomes.

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Overview

This pipeline provides tools for:
- **Metabolic gene cluster detection** using sliding window analysis
- **Homology-based searches** with BLAST against multiple databases
- **Machine learning predictions** for MGC classification
- **Phylogenetic analysis** of biosynthetic pathways
- **Integration with biological databases** (KEGG, MiBIG, PlantCyc, Ensembl)
- **Statistical analysis** and enrichment testing

## Features

### Core Functionality
- 🧬 **Genome Analysis**: Process plant genomes in FASTA format with annotation support
- 🔍 **MGC Detection**: Identify metabolic gene clusters using multiple algorithms
- 📊 **Statistical Analysis**: Enrichment testing and significance analysis
- 🤖 **Machine Learning**: Neural networks and traditional ML for cluster prediction
- 🌳 **Phylogenetic Analysis**: Evolutionary studies of biosynthetic pathways
- 📈 **Visualization**: Generate plots and reports for analysis results

### Technical Features
- ⚡ **High Performance**: Parallel processing and cluster computing support
- 🔧 **Configurable**: Flexible configuration system with environment variables
- 📝 **Comprehensive Logging**: Structured logging with performance monitoring
- 🧪 **Well Tested**: Extensive test suite with >80% coverage
- 🐳 **Containerized**: Docker support for reproducible environments
- 📚 **Well Documented**: Comprehensive API documentation and tutorials

## Installation

### Prerequisites
- Python 3.9 or higher
- BLAST+ (for sequence searches)
- HMMER (for protein domain analysis)
- Optional: SLURM (for cluster computing)

### Quick Install
```bash
# Clone the repository
git clone https://github.com/itay-mayrose/plant-mgc-analysis.git
cd plant-mgc-analysis

# Install with pip
pip install -e .

# Or install with all dependencies
pip install -e ".[dev,docs,test]"
```

### Development Install
```bash
# Clone and install in development mode
git clone https://github.com/itay-mayrose/plant-mgc-analysis.git
cd plant-mgc-analysis

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\\Scripts\\activate

# Install with all development dependencies
pip install -e ".[dev,docs,test]"

# Install pre-commit hooks
pre-commit install
```

### Docker Install
```bash
# Build Docker image
docker build -t plant-mgc-analysis .

# Run with Docker
docker run -v $(pwd)/data:/app/data plant-mgc-analysis --help
```

## Quick Start

### Basic Usage
```python
from plant_mgc_analysis import MGCAnalyzer

# Initialize analyzer
analyzer = MGCAnalyzer()

# Load genome data
genome_data = analyzer.load_genome(
    genome_file="path/to/genome.fasta",
    annotation_file="path/to/annotation.gff",
    organism_name="Arabidopsis thaliana"
)

# Configure analysis pipeline
analyzer.configure_pipeline(
    analysis_types=["sliding_window", "blast_search"],
    output_dir="./results",
    parameters={
        "sliding_window": {
            "window_size": 50000,
            "step_size": 10000,
            "min_genes": 3
        }
    }
)

# Run analysis
results = analyzer.run_pipeline()

# Get MGC candidates
candidates = analyzer.get_mgc_candidates()
print(f"Found {len(candidates)} MGC candidates")

# Generate report
report = analyzer.generate_report()
```

### Command Line Usage
```bash
# Basic MGC analysis
plant-mgc analyze --genome genome.fasta --output ./results

# With specific parameters
plant-mgc analyze \\
    --genome genome.fasta \\
    --annotation annotation.gff \\
    --organism "Arabidopsis thaliana" \\
    --analysis sliding_window blast_search \\
    --window-size 50000 \\
    --output ./results

# Run full pipeline
mgc-pipeline \\
    --input-dir ./genomes \\
    --output-dir ./results \\
    --config config.yaml \\
    --parallel-jobs 4

# Machine learning prediction
mgc-predict \\
    --candidates candidates.csv \\
    --model phytoclust \\
    --output predictions.csv
```

## Configuration

### Environment Variables
```bash
# Set environment variables
export PLANT_MGC_DATA_DIR=/path/to/data
export PLANT_MGC_OUTPUT_DIR=/path/to/output
export PLANT_MGC_CACHE_DIR=/path/to/cache

# Database settings
export PLANT_MGC_DATABASE__KEGG_API_URL=https://rest.kegg.jp
export PLANT_MGC_DATABASE__MIBIG_DATABASE_PATH=/path/to/mibig

# Compute settings
export PLANT_MGC_COMPUTE__MAX_WORKERS=8
export PLANT_MGC_COMPUTE__USE_SLURM=true
```

### Configuration File
```yaml
# config.yaml
data_dir: /path/to/data
output_dir: /path/to/output

database:
  kegg_api_url: https://rest.kegg.jp
  mibig_database_path: /path/to/mibig
  
analysis:
  default_window_size: 50000
  blast_evalue: 1e-5
  significance_threshold: 0.05
  
compute:
  max_workers: 8
  use_slurm: true
  slurm_partition: general
```

## Analysis Types

### Sliding Window Analysis
Identifies MGC candidates by scanning genomic windows for metabolic gene enrichment.

```python
results = analyzer.run_analysis(
    "sliding_window",
    parameters={
        "window_size": 50000,
        "step_size": 10000,
        "min_genes": 3,
        "significance_threshold": 0.05
    }
)
```

### BLAST Search
Performs homology searches against MGC databases.

```python
results = analyzer.run_analysis(
    "blast_search",
    parameters={
        "database": "mibig",
        "evalue": 1e-5,
        "max_targets": 100
    }
)
```

### Machine Learning
Uses trained models for MGC prediction.

```python
results = analyzer.run_analysis(
    "machine_learning",
    parameters={
        "model": "phytoclust",
        "feature_selection": True,
        "cross_validation": 5
    }
)
```

## Database Integration

### Supported Databases
- **KEGG**: Metabolic pathway information
- **MiBIG**: Known biosynthetic gene clusters
- **PlantCyc**: Plant metabolic pathways
- **Ensembl**: Plant genome annotations
- **Phytozome**: Plant genome portal

### Database Setup
```bash
# Download and setup databases
plant-mgc setup-databases --all

# Or setup specific databases
plant-mgc setup-databases --kegg --mibig
```

## Development

### Running Tests
```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=plant_mgc_analysis

# Run specific test types
pytest tests/unit/
pytest tests/integration/
```

### Code Quality
```bash
# Format code
black src/
isort src/

# Run linting
flake8 src/
mypy src/

# Run pre-commit hooks
pre-commit run --all-files
```

### Documentation
```bash
# Build documentation
cd docs/
make html

# View documentation
open docs/_build/html/index.html
```

## Project Structure

```
plant-mgc-analysis/
├── src/plant_mgc_analysis/          # Main package
│   ├── core/                        # Core functionality
│   ├── genomics/                    # Genomic analysis tools
│   ├── metabolic/                   # Metabolic pathway analysis
│   ├── machine_learning/            # ML models and predictions
│   ├── cli/                         # Command-line interface
│   ├── config/                      # Configuration management
│   ├── utils/                       # Utility functions
│   └── data/                        # Data handling
├── tests/                           # Test suite
│   ├── unit/                        # Unit tests
│   ├── integration/                 # Integration tests
│   └── fixtures/                    # Test fixtures
├── docs/                            # Documentation
├── examples/                        # Usage examples
├── scripts/                         # Utility scripts
├── data/                            # Sample data
└── docker/                          # Docker configuration
```

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

### Development Workflow
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Run the test suite
6. Submit a pull request

### Code Standards
- Follow PEP 8 style guidelines
- Add type hints to all functions
- Write comprehensive docstrings
- Maintain >80% test coverage
- Use meaningful commit messages

## Citation

If you use this software in your research, please cite:

```bibtex
@software{plant_mgc_analysis,
  title = {Plant MGC Analysis Pipeline: A Comprehensive Toolkit for Plant Biosynthetic Gene Cluster Analysis},
  author = {Mayrose, Itay and Contributors},
  year = {2024},
  url = {https://github.com/itay-mayrose/plant-mgc-analysis}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

- **Documentation**: [https://plant-mgc-analysis.readthedocs.io/](https://plant-mgc-analysis.readthedocs.io/)
- **Issues**: [GitHub Issues](https://github.com/itay-mayrose/plant-mgc-analysis/issues)
- **Discussions**: [GitHub Discussions](https://github.com/itay-mayrose/plant-mgc-analysis/discussions)

## Acknowledgments

- KEGG database for metabolic pathway information
- MiBIG database for biosynthetic gene cluster data
- PlantCyc for plant-specific metabolic pathways
- The BioPython community for sequence analysis tools
- All contributors and users of this software

## Related Projects

- [antiSMASH](https://antismash.secondarymetabolites.org/): Bacterial/fungal BGC analysis
- [ClusterFinder](https://github.com/petercr/ClusterFinder): BGC prediction tool
- [PRISM](https://github.com/magarveylab/prism-releases): Natural product discovery
- [PlantCyc](https://www.plantcyc.org/): Plant metabolic pathway database