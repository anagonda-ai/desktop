"""
PhytoClust Machine Learning Pipeline for Plant MGC Analysis.

This module provides object-oriented machine learning functionality for
metabolic gene cluster prediction using multiple ML models and feature sets.
"""

import pickle
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass, field
from pathlib import Path
import json

import numpy as np
import pandas as pd
import joblib
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import classification_report, confusion_matrix
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras.models import Sequential, load_model
    from tensorflow.keras.layers import Dense, Dropout, LSTM, Conv1D, MaxPooling1D, Flatten
    TENSORFLOW_AVAILABLE = True
except ImportError:
    TENSORFLOW_AVAILABLE = False

from ..core.base import AnalysisEngine, BioinformaticsProcessor
from ..core.types import (
    AnalysisResult,
    AnalysisType,
    GenomeData,
    GeneInfo,
    MGCCandidate,
    MLModel,
    MLFeatures,
    MLPrediction,
)
from ..core.exceptions import AnalysisError, ValidationError, ComputeError
from ..utils.file_operations import FastaProcessor, CsvProcessor
from ..utils.validation import validate_analysis_parameters


@dataclass
class FeatureSet:
    """Represents a set of features for ML training/prediction."""
    
    sequence_features: Dict[str, float] = field(default_factory=dict)
    expression_features: Dict[str, float] = field(default_factory=dict)
    domain_features: Dict[str, float] = field(default_factory=dict)
    genomic_features: Dict[str, float] = field(default_factory=dict)
    metabolic_features: Dict[str, float] = field(default_factory=dict)
    
    def to_vector(self, feature_names: List[str]) -> np.ndarray:
        """Convert to feature vector using specified feature names."""
        all_features = {
            **self.sequence_features,
            **self.expression_features,
            **self.domain_features,
            **self.genomic_features,
            **self.metabolic_features,
        }
        
        vector = []
        for name in feature_names:
            vector.append(all_features.get(name, 0.0))
        
        return np.array(vector)
    
    def get_all_features(self) -> Dict[str, float]:
        """Get all features as a single dictionary."""
        return {
            **self.sequence_features,
            **self.expression_features,
            **self.domain_features,
            **self.genomic_features,
            **self.metabolic_features,
        }


@dataclass
class MLModelResult:
    """Results from ML model prediction."""
    
    prediction: Union[int, float, str]
    confidence: float
    probabilities: Optional[np.ndarray] = None
    features_used: List[str] = field(default_factory=list)
    model_name: str = ""
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        result = {
            "prediction": self.prediction,
            "confidence": self.confidence,
            "model_name": self.model_name,
            "features_used": self.features_used,
        }
        
        if self.probabilities is not None:
            result["probabilities"] = self.probabilities.tolist()
        
        return result


class FeatureExtractor(BioinformaticsProcessor):
    """Extracts features from biological sequences and annotations."""
    
    def __init__(self, **kwargs):
        """Initialize feature extractor."""
        super().__init__(**kwargs)
        
        # Load expression data if available
        self.expression_data: Optional[pd.DataFrame] = None
        self.hmmer_profiles: Optional[Path] = None
        
        # Core gene patterns for metabolic pathways
        self.core_gene_patterns = [
            "CYP", "PKS", "TPS", "NRPS", "UGT", "ABC", "MT", "AT", "KS"
        ]
    
    def load_expression_data(self, expression_file: Path) -> None:
        """Load gene expression data."""
        processor = CsvProcessor()
        self.expression_data = processor.read_file(expression_file)
        self.expression_data.set_index(self.expression_data.columns[0], inplace=True)
        
        self.logger.info(f"Loaded expression data for {len(self.expression_data)} genes")
    
    def load_hmmer_profiles(self, profile_path: Path) -> None:
        """Load HMMER profile database path."""
        self.hmmer_profiles = profile_path
        self.logger.info(f"Set HMMER profiles path: {profile_path}")
    
    def extract_sequence_features(self, sequence: str, gene_id: str = "") -> Dict[str, float]:
        """
        Extract sequence-based features.
        
        Args:
            sequence: Protein or nucleotide sequence
            gene_id: Gene identifier
            
        Returns:
            Dictionary of sequence features
        """
        if not sequence:
            return {}
        
        seq = sequence.upper()
        length = len(seq)
        
        features = {
            "sequence_length": float(length),
            "gc_content": (seq.count("G") + seq.count("C")) / length if length > 0 else 0.0,
        }
        
        # Amino acid composition (if protein sequence)
        if all(c in "ACDEFGHIKLMNPQRSTVWY" for c in seq[:100]):  # Check if protein
            aa_counts = {aa: seq.count(aa) / length for aa in "ACDEFGHIKLMNPQRSTVWY"}
            features.update({f"aa_{aa}": count for aa, count in aa_counts.items()})
            
            # Secondary structure indicators
            hydrophobic = sum(seq.count(aa) for aa in "AILMFPWYV") / length
            charged = sum(seq.count(aa) for aa in "DEKR") / length
            polar = sum(seq.count(aa) for aa in "STNQHCY") / length
            
            features.update({
                "hydrophobic_fraction": hydrophobic,
                "charged_fraction": charged,
                "polar_fraction": polar,
            })
        
        # Core gene indicator
        features["is_core_gene"] = float(any(
            pattern in gene_id.upper() for pattern in self.core_gene_patterns
        ))
        
        return features
    
    def extract_expression_features(self, gene_id: str) -> Dict[str, float]:
        """Extract expression-based features."""
        features = {}
        
        if self.expression_data is not None and gene_id in self.expression_data.index:
            expression_values = self.expression_data.loc[gene_id].values
            
            features.update({
                "mean_expression": float(np.mean(expression_values)),
                "max_expression": float(np.max(expression_values)),
                "min_expression": float(np.min(expression_values)),
                "std_expression": float(np.std(expression_values)),
                "cv_expression": float(np.std(expression_values) / np.mean(expression_values))
                    if np.mean(expression_values) > 0 else 0.0,
            })
        else:
            # Default values if no expression data
            features.update({
                "mean_expression": 0.0,
                "max_expression": 0.0,
                "min_expression": 0.0,
                "std_expression": 0.0,
                "cv_expression": 0.0,
            })
        
        return features
    
    def extract_domain_features(self, sequence: str, gene_id: str = "") -> Dict[str, float]:
        """
        Extract protein domain features using HMMER.
        
        Args:
            sequence: Protein sequence
            gene_id: Gene identifier
            
        Returns:
            Dictionary of domain features
        """
        features = {}
        
        if not self.hmmer_profiles or not sequence:
            return features
        
        try:
            # Create temporary files
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
                temp_fasta.write(f">{gene_id}\n{sequence}\n")
                temp_fasta_path = temp_fasta.name
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.out', delete=False) as temp_out:
                temp_out_path = temp_out.name
            
            # Run HMMER
            cmd = [
                "hmmscan",
                "--cpu", "1",
                "--domtblout", temp_out_path,
                str(self.hmmer_profiles),
                temp_fasta_path
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                # Parse HMMER results
                domain_scores = []
                with open(temp_out_path, 'r') as f:
                    for line in f:
                        if not line.startswith('#'):
                            parts = line.split()
                            if len(parts) >= 14:
                                score = float(parts[13])  # Domain score
                                domain_scores.append(score)
                
                if domain_scores:
                    features.update({
                        "domain_count": float(len(domain_scores)),
                        "max_domain_score": float(max(domain_scores)),
                        "mean_domain_score": float(np.mean(domain_scores)),
                        "total_domain_score": float(sum(domain_scores)),
                    })
                else:
                    features.update({
                        "domain_count": 0.0,
                        "max_domain_score": 0.0,
                        "mean_domain_score": 0.0,
                        "total_domain_score": 0.0,
                    })
            
            # Cleanup
            Path(temp_fasta_path).unlink()
            Path(temp_out_path).unlink()
            
        except Exception as e:
            self.logger.warning(f"Failed to extract domain features for {gene_id}: {e}")
            features.update({
                "domain_count": 0.0,
                "max_domain_score": 0.0,
                "mean_domain_score": 0.0,
                "total_domain_score": 0.0,
            })
        
        return features
    
    def extract_genomic_features(self, gene: GeneInfo, neighbors: List[GeneInfo] = None) -> Dict[str, float]:
        """
        Extract genomic context features.
        
        Args:
            gene: Gene information
            neighbors: Neighboring genes
            
        Returns:
            Dictionary of genomic features
        """
        features = {
            "gene_length": float(gene.length),
            "strand": 1.0 if gene.strand == "+" else -1.0,
        }
        
        if neighbors:
            # Distance to nearest neighbors
            distances = []
            for neighbor in neighbors:
                if neighbor.gene_id != gene.gene_id:
                    distance = min(
                        abs(gene.start - neighbor.end),
                        abs(gene.end - neighbor.start)
                    )
                    distances.append(distance)
            
            if distances:
                features.update({
                    "min_neighbor_distance": float(min(distances)),
                    "mean_neighbor_distance": float(np.mean(distances)),
                    "neighbor_count": float(len(neighbors) - 1),
                })
            else:
                features.update({
                    "min_neighbor_distance": 0.0,
                    "mean_neighbor_distance": 0.0,
                    "neighbor_count": 0.0,
                })
        
        return features
    
    def extract_metabolic_features(self, gene: GeneInfo) -> Dict[str, float]:
        """Extract metabolic pathway features."""
        features = {
            "has_ec_number": 1.0 if hasattr(gene, 'ec_number') and gene.ec_number else 0.0,
            "has_pathway": 1.0 if gene.pathway_ids else 0.0,
            "pathway_count": float(len(gene.pathway_ids)),
        }
        
        # EC number analysis
        if hasattr(gene, 'ec_number') and gene.ec_number:
            ec_parts = gene.ec_number.split('.')
            if len(ec_parts) >= 2:
                features.update({
                    "ec_class": float(ec_parts[0]) if ec_parts[0].isdigit() else 0.0,
                    "ec_subclass": float(ec_parts[1]) if ec_parts[1].isdigit() else 0.0,
                })
        
        return features
    
    def extract_features(self, gene: GeneInfo, sequence: str = "", neighbors: List[GeneInfo] = None) -> FeatureSet:
        """
        Extract complete feature set for a gene.
        
        Args:
            gene: Gene information
            sequence: Protein sequence
            neighbors: Neighboring genes
            
        Returns:
            Complete feature set
        """
        return FeatureSet(
            sequence_features=self.extract_sequence_features(sequence, gene.gene_id),
            expression_features=self.extract_expression_features(gene.gene_id),
            domain_features=self.extract_domain_features(sequence, gene.gene_id),
            genomic_features=self.extract_genomic_features(gene, neighbors),
            metabolic_features=self.extract_metabolic_features(gene),
        )
    
    def validate_input(self, data: Any) -> None:
        """Validate input data."""
        if not isinstance(data, (GeneInfo, list, str)):
            raise ValidationError(
                "Input must be GeneInfo object, gene list, or sequence",
                field_name="input_data",
                field_value=type(data).__name__
            )
    
    def process(self, data: Any, **kwargs) -> FeatureSet:
        """Process feature extraction."""
        self.validate_input(data)
        
        if isinstance(data, GeneInfo):
            sequence = kwargs.get("sequence", "")
            neighbors = kwargs.get("neighbors", [])
            return self.extract_features(data, sequence, neighbors)
        
        elif isinstance(data, str):
            # Direct sequence processing
            gene_id = kwargs.get("gene_id", "unknown")
            return FeatureSet(
                sequence_features=self.extract_sequence_features(data, gene_id),
                expression_features=self.extract_expression_features(gene_id),
            )
        
        return FeatureSet()


class MLModelManager:
    """Manages multiple ML models for different prediction tasks."""
    
    def __init__(self):
        """Initialize model manager."""
        self.models: Dict[str, Any] = {}
        self.scalers: Dict[str, StandardScaler] = {}
        self.feature_names: Dict[str, List[str]] = {}
        self.label_encoders: Dict[str, LabelEncoder] = {}
    
    def add_model(
        self,
        name: str,
        model: Any,
        scaler: Optional[StandardScaler] = None,
        feature_names: Optional[List[str]] = None,
        label_encoder: Optional[LabelEncoder] = None,
    ) -> None:
        """Add a model to the manager."""
        self.models[name] = model
        
        if scaler:
            self.scalers[name] = scaler
        if feature_names:
            self.feature_names[name] = feature_names
        if label_encoder:
            self.label_encoders[name] = label_encoder
    
    def load_model(self, name: str, model_path: Path, model_type: str = "sklearn") -> None:
        """
        Load a model from file.
        
        Args:
            name: Model name
            model_path: Path to model file
            model_type: Type of model (sklearn, tensorflow, custom)
        """
        try:
            if model_type == "sklearn":
                model = joblib.load(model_path)
            elif model_type == "tensorflow" and TENSORFLOW_AVAILABLE:
                model = load_model(model_path)
            elif model_type == "pickle":
                with open(model_path, 'rb') as f:
                    model = pickle.load(f)
            else:
                raise ValueError(f"Unsupported model type: {model_type}")
            
            self.models[name] = model
            
        except Exception as e:
            raise ComputeError(
                f"Failed to load model {name}: {e}",
                computation_type="model_loading",
                resource_type="model_file"
            ) from e
    
    def save_model(self, name: str, model_path: Path, model_type: str = "sklearn") -> None:
        """Save a model to file."""
        if name not in self.models:
            raise ValueError(f"Model {name} not found")
        
        model = self.models[name]
        
        try:
            if model_type == "sklearn":
                joblib.dump(model, model_path)
            elif model_type == "tensorflow" and TENSORFLOW_AVAILABLE:
                model.save(model_path)
            elif model_type == "pickle":
                with open(model_path, 'wb') as f:
                    pickle.dump(model, f)
            
        except Exception as e:
            raise ComputeError(
                f"Failed to save model {name}: {e}",
                computation_type="model_saving",
                resource_type="model_file"
            ) from e
    
    def predict(self, name: str, features: Union[FeatureSet, np.ndarray]) -> MLModelResult:
        """Make prediction using specified model."""
        if name not in self.models:
            raise ValueError(f"Model {name} not found")
        
        model = self.models[name]
        
        # Prepare features
        if isinstance(features, FeatureSet):
            if name in self.feature_names:
                feature_vector = features.to_vector(self.feature_names[name])
            else:
                feature_vector = np.array(list(features.get_all_features().values()))
        else:
            feature_vector = features
        
        # Reshape for single sample
        if feature_vector.ndim == 1:
            feature_vector = feature_vector.reshape(1, -1)
        
        # Scale features if scaler available
        if name in self.scalers:
            feature_vector = self.scalers[name].transform(feature_vector)
        
        # Make prediction
        try:
            if hasattr(model, 'predict_proba'):
                probabilities = model.predict_proba(feature_vector)[0]
                prediction = model.predict(feature_vector)[0]
                confidence = float(np.max(probabilities))
            else:
                prediction = model.predict(feature_vector)[0]
                probabilities = None
                confidence = 1.0  # Default confidence for non-probabilistic models
            
            # Decode labels if encoder available
            if name in self.label_encoders:
                prediction = self.label_encoders[name].inverse_transform([prediction])[0]
            
            return MLModelResult(
                prediction=prediction,
                confidence=confidence,
                probabilities=probabilities,
                features_used=self.feature_names.get(name, []),
                model_name=name,
            )
            
        except Exception as e:
            raise ComputeError(
                f"Prediction failed for model {name}: {e}",
                computation_type="model_prediction",
                resource_type="model"
            ) from e


class PhytoClustAnalyzer(AnalysisEngine):
    """
    PhytoClust machine learning analyzer for MGC prediction.
    
    This class implements a comprehensive ML pipeline for metabolic gene
    cluster detection using multiple models and feature types.
    """
    
    def __init__(self, **kwargs):
        """Initialize PhytoClust analyzer."""
        super().__init__(
            analysis_type=AnalysisType.MACHINE_LEARNING,
            **kwargs
        )
        
        self.feature_extractor = FeatureExtractor(**kwargs)
        self.model_manager = MLModelManager()
        
        # Default parameters
        self.cluster_max_distance = 10000  # Maximum distance between genes in cluster
        self.min_genes_per_cluster = 2     # Minimum genes required for cluster
        self.confidence_threshold = 0.7    # Minimum confidence for positive prediction
    
    def validate_input(self, genome_data: GenomeData) -> None:
        """Validate input genome data."""
        if not genome_data.protein_file:
            raise ValidationError(
                "Protein file required for ML analysis",
                field_name="protein_file",
                field_value=None
            )
    
    def load_models(self, model_config: Dict[str, Dict[str, Any]]) -> None:
        """
        Load ML models from configuration.
        
        Args:
            model_config: Dictionary with model configurations
        """
        for model_name, config in model_config.items():
            model_path = Path(config["path"])
            model_type = config.get("type", "sklearn")
            
            self.model_manager.load_model(model_name, model_path, model_type)
            
            # Load associated files
            if "scaler_path" in config:
                scaler = joblib.load(config["scaler_path"])
                self.model_manager.scalers[model_name] = scaler
            
            if "features_path" in config:
                with open(config["features_path"], 'r') as f:
                    features = json.load(f)
                self.model_manager.feature_names[model_name] = features
        
        self.logger.info(f"Loaded {len(model_config)} ML models")
    
    def predict_gene_clusters(self, genes: List[GeneInfo], sequences: Dict[str, str]) -> List[MGCCandidate]:
        """
        Predict metabolic gene clusters using ML models.
        
        Args:
            genes: List of genes to analyze
            sequences: Dictionary mapping gene IDs to protein sequences
            
        Returns:
            List of predicted MGC candidates
        """
        candidates = []
        
        # Extract features for all genes
        gene_features = {}
        metabolic_scores = {}
        
        for gene in genes:
            sequence = sequences.get(gene.gene_id, "")
            
            # Get neighboring genes (within 50kb)
            neighbors = [
                g for g in genes
                if g.chromosome == gene.chromosome and
                abs(g.start - gene.start) <= 50000
            ]
            
            # Extract features
            features = self.feature_extractor.extract_features(gene, sequence, neighbors)
            gene_features[gene.gene_id] = features
            
            # Predict if gene is metabolic
            if "metabolic_classifier" in self.model_manager.models:
                result = self.model_manager.predict("metabolic_classifier", features)
                metabolic_scores[gene.gene_id] = result.confidence if result.prediction else 0.0
        
        # Group genes by chromosome
        genes_by_chr = {}
        for gene in genes:
            if gene.chromosome not in genes_by_chr:
                genes_by_chr[gene.chromosome] = []
            genes_by_chr[gene.chromosome].append(gene)
        
        # Detect clusters on each chromosome
        cluster_id = 0
        for chromosome, chr_genes in genes_by_chr.items():
            # Filter metabolic genes
            metabolic_genes = [
                gene for gene in chr_genes
                if metabolic_scores.get(gene.gene_id, 0.0) > self.confidence_threshold
            ]
            
            if len(metabolic_genes) < self.min_genes_per_cluster:
                continue
            
            # Sort by position
            metabolic_genes.sort(key=lambda g: g.start)
            
            # Find clusters
            current_cluster = []
            
            for gene in metabolic_genes:
                if not current_cluster:
                    current_cluster.append(gene)
                else:
                    # Check distance to last gene in cluster
                    last_gene = current_cluster[-1]
                    distance = gene.start - last_gene.end
                    
                    if distance <= self.cluster_max_distance:
                        current_cluster.append(gene)
                    else:
                        # End current cluster and start new one
                        if len(current_cluster) >= self.min_genes_per_cluster:
                            candidate = self._create_candidate_from_cluster(
                                current_cluster, cluster_id, gene_features, metabolic_scores
                            )
                            candidates.append(candidate)
                            cluster_id += 1
                        
                        current_cluster = [gene]
            
            # Add final cluster
            if len(current_cluster) >= self.min_genes_per_cluster:
                candidate = self._create_candidate_from_cluster(
                    current_cluster, cluster_id, gene_features, metabolic_scores
                )
                candidates.append(candidate)
                cluster_id += 1
        
        return candidates
    
    def _create_candidate_from_cluster(
        self,
        cluster_genes: List[GeneInfo],
        cluster_id: int,
        gene_features: Dict[str, FeatureSet],
        metabolic_scores: Dict[str, float],
    ) -> MGCCandidate:
        """Create MGC candidate from gene cluster."""
        
        # Calculate cluster boundaries
        cluster_start = min(gene.start for gene in cluster_genes)
        cluster_end = max(gene.end for gene in cluster_genes)
        
        # Calculate confidence score
        scores = [metabolic_scores.get(gene.gene_id, 0.0) for gene in cluster_genes]
        confidence_score = np.mean(scores)
        
        # Predict pathway if pathway model available
        pathway_prediction = None
        if "pathway_classifier" in self.model_manager.models:
            # Aggregate features for cluster
            cluster_features = self._aggregate_cluster_features(cluster_genes, gene_features)
            pathway_result = self.model_manager.predict("pathway_classifier", cluster_features)
            pathway_prediction = pathway_result.prediction
        
        return MGCCandidate(
            cluster_id=f"ML_{cluster_id:04d}",
            chromosome=cluster_genes[0].chromosome,
            start=cluster_start,
            end=cluster_end,
            genes=cluster_genes,
            organism="unknown",  # Will be set by caller
            analysis_method="phytoclust_ml",
            confidence_score=confidence_score,
            pathway_name=pathway_prediction,
            validation_data={
                "metabolic_scores": {gene.gene_id: metabolic_scores.get(gene.gene_id, 0.0) 
                                   for gene in cluster_genes},
                "cluster_length": cluster_end - cluster_start,
                "gene_count": len(cluster_genes),
            }
        )
    
    def _aggregate_cluster_features(
        self,
        cluster_genes: List[GeneInfo],
        gene_features: Dict[str, FeatureSet],
    ) -> FeatureSet:
        """Aggregate features across genes in a cluster."""
        
        all_features = {}
        feature_counts = {}
        
        for gene in cluster_genes:
            if gene.gene_id in gene_features:
                features_dict = gene_features[gene.gene_id].get_all_features()
                
                for feature_name, value in features_dict.items():
                    if feature_name not in all_features:
                        all_features[feature_name] = 0.0
                        feature_counts[feature_name] = 0
                    
                    all_features[feature_name] += value
                    feature_counts[feature_name] += 1
        
        # Calculate means
        aggregated_features = {}
        for feature_name, total_value in all_features.items():
            count = feature_counts[feature_name]
            aggregated_features[feature_name] = total_value / count if count > 0 else 0.0
        
        # Add cluster-specific features
        aggregated_features.update({
            "cluster_gene_count": float(len(cluster_genes)),
            "cluster_length": float(max(gene.end for gene in cluster_genes) - 
                                   min(gene.start for gene in cluster_genes)),
        })
        
        return FeatureSet(metabolic_features=aggregated_features)
    
    def analyze(
        self,
        genome_data: GenomeData,
        parameters: Optional[Dict[str, Any]] = None,
    ) -> AnalysisResult:
        """
        Perform ML-based MGC analysis.
        
        Args:
            genome_data: Genome data to analyze
            parameters: Optional analysis parameters
            
        Returns:
            Analysis result
        """
        if parameters:
            self.validate_parameters(parameters)
        
        # Load models if specified
        if "model_config" in (parameters or {}):
            self.load_models(parameters["model_config"])
        
        # Load protein sequences
        processor = FastaProcessor()
        sequences = processor.read_file(genome_data.protein_file)
        sequence_dict = {seq.id: str(seq.seq) for seq in sequences}
        
        # Load gene annotations
        if genome_data.annotation_file:
            from ..utils.file_operations import GffProcessor
            gff_processor = GffProcessor()
            genes = gff_processor.read_file(genome_data.annotation_file)
        else:
            # Create gene objects from sequences
            genes = [
                GeneInfo(
                    gene_id=seq.id,
                    chromosome="unknown",
                    start=0,
                    end=len(seq.seq),
                    strand="+",
                )
                for seq in sequences
            ]
        
        # Predict clusters
        candidates = self.predict_gene_clusters(genes, sequence_dict)
        
        # Set organism name
        for candidate in candidates:
            candidate.organism = genome_data.organism_name
        
        # Create analysis result
        result = AnalysisResult(
            analysis_id=f"{self.session_id}_phytoclust",
            analysis_type=AnalysisType.MACHINE_LEARNING,
            input_data={
                "genome_file": str(genome_data.genome_file),
                "protein_file": str(genome_data.protein_file),
                "organism": genome_data.organism_name,
            },
            results={
                "mgc_candidates": candidates,
                "total_genes": len(genes),
                "candidate_count": len(candidates),
                "models_used": list(self.model_manager.models.keys()),
            },
            parameters=parameters or {},
            status="success",
        )
        
        self.logger.info(
            f"PhytoClust analysis complete: {len(candidates)} candidates predicted"
        )
        
        return result
    
    def process(self, data: GenomeData, **kwargs) -> AnalysisResult:
        """Process genome data with PhytoClust analysis."""
        return self.analyze(data, kwargs.get('parameters'))
    
    def validate_parameters(self, parameters: Dict[str, Any]) -> None:
        """Validate analysis parameters."""
        super().validate_parameters(parameters)
        
        validated = validate_analysis_parameters("machine_learning", parameters)
        
        # Update instance parameters
        self.cluster_max_distance = validated.get("cluster_max_distance", self.cluster_max_distance)
        self.min_genes_per_cluster = validated.get("min_genes_per_cluster", self.min_genes_per_cluster)
        self.confidence_threshold = validated.get("confidence_threshold", self.confidence_threshold)