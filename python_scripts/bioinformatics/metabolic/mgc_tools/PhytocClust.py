import os
import json
import numpy as np
import pandas as pd
from flask import Flask, request, jsonify
import tempfile
import requests
import subprocess
from Bio import SeqIO
from Bio.SearchIO import parse as hmm_parse
from sklearn.preprocessing import StandardScaler
import joblib

# -------------------------------
# Object-Oriented AI-Optimized Gene Cluster Analyzer with HMMER Annotation & Pathway Prediction
# -------------------------------
class GeneClusterAnalyzer:
    def __init__(self, model_path, pathway_model_path, cluster_max_distance=10000, min_genes_per_cluster=3):
        self.model = joblib.load(model_path)
        self.pathway_model = joblib.load(pathway_model_path)
        self.scaler = StandardScaler()
        self.cluster_max_distance = cluster_max_distance
        self.min_genes_per_cluster = min_genes_per_cluster

    def save_uploaded_file(self, uploaded_file, output_dir):
        file_path = os.path.join(output_dir, uploaded_file.filename)
        uploaded_file.save(file_path)
        return file_path

    def parse_fasta(self, fasta_file):
        return list(SeqIO.parse(fasta_file, "fasta"))

    def predict_genes(self, fasta_file):
        records = self.parse_fasta(fasta_file)
        predicted_genes = []
        for record in records:
            for strand, start in [(+1, 0), (-1, 1)]:
                for i in range(start, len(record.seq), 3):
                    codon = record.seq[i:i+3]
                    if codon in ["ATG"]:
                        for j in range(i, len(record.seq)-3, 3):
                            stop_codon = record.seq[j:j+3]
                            if stop_codon in ["TAA", "TAG", "TGA"]:
                                gene_seq = str(record.seq[i:j+3])
                                predicted_genes.append({
                                    "id": f"{record.id}_gene_{i}",
                                    "sequence": gene_seq,
                                    "start": i,
                                    "end": j+3,
                                    "strand": strand,
                                    "length": j+3 - i
                                })
                                break
        return predicted_genes

    def annotate_genes_with_hmmer(self, genes, hmm_db):
        """Annotates genes using HMMER model."""
        annotated_genes = []
        for gene in genes:
            tmp_fasta = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta")
            tmp_fasta.write(f">{gene['id']}\n{gene['sequence']}\n".encode())
            tmp_fasta.close()
            
            cmd = f"hmmscan --tblout {tmp_fasta.name}.out {hmm_db} {tmp_fasta.name}"
            subprocess.run(cmd, shell=True)
            
            with open(f"{tmp_fasta.name}.out", "r") as result_file:
                for line in result_file:
                    if not line.startswith("#"):
                        parts = line.split()
                        if len(parts) > 0:
                            gene['annotation'] = parts[0]
                            break
            
            annotated_genes.append(gene)
        return annotated_genes
    
    def predict_pathways(self, annotated_genes):
        """Predicts metabolic pathways based on AI model."""
        pathway_predictions = {}
        for gene in annotated_genes:
            features = np.array([[len(gene['sequence'])]]).reshape(1, -1)
            predicted_pathway = self.pathway_model.predict(features)[0]
            pathway_predictions[gene['id']] = predicted_pathway
            gene['pathway'] = predicted_pathway
        return annotated_genes

    def detect_metabolic_gene_clusters(self, annotated_genes):
        """Groups genes into clusters based on shared pathways and spatial proximity."""
        clusters = []
        sorted_genes = sorted(annotated_genes, key=lambda x: x['start'])
        current_cluster = []
        current_pathway = None
        for i, gene in enumerate(sorted_genes):
            if not current_cluster:
                current_cluster.append(gene)
                current_pathway = gene['pathway']
            else:
                last_gene = current_cluster[-1]
                if gene['start'] - last_gene['end'] <= self.cluster_max_distance and gene['pathway'] == current_pathway:
                    current_cluster.append(gene)
                else:
                    if len(current_cluster) >= self.min_genes_per_cluster:
                        clusters.append(current_cluster)
                    current_cluster = [gene]
                    current_pathway = gene['pathway']
        if len(current_cluster) >= self.min_genes_per_cluster:
            clusters.append(current_cluster)
        return clusters

# -------------------------------
# Flask API for User Interaction
# -------------------------------
app = Flask(__name__)
analyzer = GeneClusterAnalyzer("phytoClust_optimized_model.pkl", "pathway_prediction_model.pkl", cluster_max_distance=10000, min_genes_per_cluster=3)

@app.route('/predict', methods=['POST'])
def predict_clusters():
    genome_file = request.files['genome']
    hmm_db = request.form['hmm_db']
    
    with tempfile.TemporaryDirectory() as output_dir:
        genome_path = analyzer.save_uploaded_file(genome_file, output_dir)
        predicted_genes = analyzer.predict_genes(genome_path)
        annotated_genes = analyzer.annotate_genes_with_hmmer(predicted_genes, hmm_db)
        pathway_assigned_genes = analyzer.predict_pathways(annotated_genes)
        clusters = analyzer.detect_metabolic_gene_clusters(pathway_assigned_genes)
    
    response = {
        "clusters": clusters
    }
    return jsonify(response)

if __name__ == '__main__':
    app.run(debug=True)
