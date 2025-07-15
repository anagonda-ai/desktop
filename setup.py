"""
Setup script for Plant MGC Analysis Pipeline.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text() if (this_directory / "README.md").exists() else ""

setup(
    name="plant-mgc-analysis",
    version="1.0.0",
    description="Unified Plant Metabolic Gene Cluster Analysis Pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Plant MGC Analysis Team",
    author_email="user@example.com",
    url="https://github.com/your-org/plant-mgc-analysis",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.8",
    install_requires=[
        "click>=8.0.0",
        "loguru>=0.6.0",
        "pydantic>=2.0.0",
        "pydantic-settings>=2.0.0",
        "pandas>=1.5.0",
        "numpy>=1.21.0",
        "biopython>=1.79",
        "requests>=2.28.0",
        "scipy>=1.9.0",
        "scikit-learn>=1.1.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.11.0",
        "tqdm>=4.64.0",
        "psutil>=5.9.0",
        "aiofiles>=0.8.0",
        "httpx>=0.23.0",
        "tabulate>=0.9.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "pytest-asyncio>=0.20.0",
            "black>=22.0.0",
            "isort>=5.10.0",
            "flake8>=5.0.0",
            "mypy>=0.991",
            "pre-commit>=2.20.0",
        ],
        "docs": [
            "mkdocs>=1.4.0",
            "mkdocs-material>=8.5.0",
            "mkdocstrings>=0.19.0",
        ],
        "cluster": [
            "dask[complete]>=2022.10.0",
            "distributed>=2022.10.0",
            "slurm-jobqueue>=0.8.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "plant-mgc-analysis=plant_mgc_analysis.cli_main:cli",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="bioinformatics genomics metabolic-pathways plant-biology",
    project_urls={
        "Bug Reports": "https://github.com/your-org/plant-mgc-analysis/issues",
        "Source": "https://github.com/your-org/plant-mgc-analysis",
        "Documentation": "https://plant-mgc-analysis.readthedocs.io/",
    },
)