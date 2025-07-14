"""
Test configuration and fixtures for Plant MGC Analysis Pipeline.

This module provides common test fixtures and configuration for the test suite.
"""

import tempfile
from pathlib import Path
from typing import Generator

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from plant_mgc_analysis import MGCAnalyzer
from plant_mgc_analysis.core.types import GenomeData, GeneInfo, MGCCandidate
from plant_mgc_analysis.config.settings import Settings


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """Provide a temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


@pytest.fixture
def sample_fasta_file(temp_dir: Path) -> Path:
    """Create a sample FASTA file for testing."""
    sequences = [
        SeqRecord(
            Seq("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"),
            id="chr1",
            description="Chromosome 1"
        ),
        SeqRecord(
            Seq("GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
            id="chr2", 
            description="Chromosome 2"
        ),
    ]
    
    fasta_file = temp_dir / "test_genome.fasta"
    with open(fasta_file, "w") as f:
        SeqIO.write(sequences, f, "fasta")
    
    return fasta_file


@pytest.fixture
def sample_gff_file(temp_dir: Path) -> Path:
    """Create a sample GFF file for testing."""
    gff_content = """##gff-version 3
chr1	test	gene	100	500	.	+	.	ID=gene1;Name=test_gene1
chr1	test	mRNA	100	500	.	+	.	ID=mRNA1;Parent=gene1
chr1	test	exon	100	200	.	+	.	ID=exon1;Parent=mRNA1
chr1	test	exon	300	500	.	+	.	ID=exon2;Parent=mRNA1
chr2	test	gene	1000	1500	.	-	.	ID=gene2;Name=test_gene2
chr2	test	mRNA	1000	1500	.	-	.	ID=mRNA2;Parent=gene2
chr2	test	exon	1000	1500	.	-	.	ID=exon3;Parent=mRNA2
"""
    
    gff_file = temp_dir / "test_annotation.gff"
    with open(gff_file, "w") as f:
        f.write(gff_content)
    
    return gff_file


@pytest.fixture
def sample_protein_file(temp_dir: Path) -> Path:
    """Create a sample protein FASTA file for testing."""
    sequences = [
        SeqRecord(
            Seq("MKLAVLVLLVVTASTPEAARKFLKDQVDLVGVYGTEGQSSTISNYGSGD"),
            id="protein1",
            description="Test protein 1"
        ),
        SeqRecord(
            Seq("MKLVFLVLLVVTASTPEAARKFLKDQVDLVGVYGTEGQSSTISNYGSGD"),
            id="protein2",
            description="Test protein 2"
        ),
    ]
    
    protein_file = temp_dir / "test_proteins.fasta"
    with open(protein_file, "w") as f:
        SeqIO.write(sequences, f, "fasta")
    
    return protein_file


@pytest.fixture
def sample_genome_data(
    sample_fasta_file: Path,
    sample_gff_file: Path,
    sample_protein_file: Path,
) -> GenomeData:
    """Create sample genome data for testing."""
    return GenomeData(
        organism_name="Test organism",
        genome_file=sample_fasta_file,
        annotation_file=sample_gff_file,
        protein_file=sample_protein_file,
        species="Test species",
        strain="Test strain",
        metadata={"test": "data"},
    )


@pytest.fixture
def sample_gene_info() -> GeneInfo:
    """Create sample gene info for testing."""
    return GeneInfo(
        gene_id="test_gene_1",
        chromosome="chr1",
        start=100,
        end=500,
        strand="+",
        gene_name="test_gene",
        protein_id="test_protein_1",
        description="Test gene description",
        functional_annotation="Test function",
        pathway_ids=["pathway1", "pathway2"],
        metabolic_role="biosynthesis",
    )


@pytest.fixture
def sample_mgc_candidate(sample_gene_info: GeneInfo) -> MGCCandidate:
    """Create sample MGC candidate for testing."""
    return MGCCandidate(
        cluster_id="test_cluster_1",
        chromosome="chr1",
        start=100,
        end=2000,
        genes=[sample_gene_info],
        organism="Test organism",
        analysis_method="sliding_window",
        confidence_score=0.85,
        pathway_id="test_pathway",
        pathway_name="Test pathway",
        metabolic_class="terpenoid",
        validation_data={"test": "validation"},
    )


@pytest.fixture
def test_settings(temp_dir: Path) -> Settings:
    """Create test settings configuration."""
    return Settings(
        data_dir=temp_dir / "data",
        output_dir=temp_dir / "output",
        cache_dir=temp_dir / "cache",
        debug=True,
    )


@pytest.fixture
def mgc_analyzer(test_settings: Settings) -> MGCAnalyzer:
    """Create MGC analyzer instance for testing."""
    return MGCAnalyzer(config=test_settings.model_dump())


@pytest.fixture
def mock_blast_results():
    """Create mock BLAST results for testing."""
    return [
        {
            "query_id": "query1",
            "subject_id": "subject1",
            "identity": 95.0,
            "evalue": 1e-10,
            "bit_score": 150.0,
        },
        {
            "query_id": "query2",
            "subject_id": "subject2",
            "identity": 80.0,
            "evalue": 1e-5,
            "bit_score": 100.0,
        },
    ]


@pytest.fixture
def mock_kegg_response():
    """Create mock KEGG API response for testing."""
    return {
        "pathway_id": "ko00900",
        "pathway_name": "Terpenoid backbone biosynthesis",
        "genes": [
            {
                "gene_id": "K01662",
                "gene_name": "HMG-CoA reductase",
                "ec_number": "1.1.1.34",
            },
            {
                "gene_id": "K00626",
                "gene_name": "acetyl-CoA C-acetyltransferase",
                "ec_number": "2.3.1.9",
            },
        ],
    }


# Test markers for different test categories
def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line("markers", "unit: mark test as unit test")
    config.addinivalue_line("markers", "integration: mark test as integration test")
    config.addinivalue_line("markers", "slow: mark test as slow running")
    config.addinivalue_line("markers", "network: mark test as requiring network access")
    config.addinivalue_line("markers", "gpu: mark test as requiring GPU")


# Auto-use fixtures for specific test types
@pytest.fixture(autouse=True)
def setup_test_environment(request, temp_dir: Path):
    """Setup test environment for each test."""
    # Create necessary directories
    (temp_dir / "data").mkdir(exist_ok=True)
    (temp_dir / "output").mkdir(exist_ok=True)
    (temp_dir / "cache").mkdir(exist_ok=True)
    
    # Set test-specific environment variables if needed
    import os
    os.environ["PLANT_MGC_DATA_DIR"] = str(temp_dir / "data")
    os.environ["PLANT_MGC_OUTPUT_DIR"] = str(temp_dir / "output")
    os.environ["PLANT_MGC_CACHE_DIR"] = str(temp_dir / "cache")
    
    yield
    
    # Cleanup after test
    for key in ["PLANT_MGC_DATA_DIR", "PLANT_MGC_OUTPUT_DIR", "PLANT_MGC_CACHE_DIR"]:
        os.environ.pop(key, None)