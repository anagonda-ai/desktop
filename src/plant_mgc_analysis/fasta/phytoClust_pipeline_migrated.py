from ..config.settings import get_settings
from ..core.base import BioinformaticsProcessor
from ..core.exceptions import ProcessingError, ValidationError
from ..utils.file_operations import file_manager
from loguru import logger
import os
import json
import subprocess
import numpy as np
import pandas as pd
import joblib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tensorflow.keras.models import load_model
from sklearn.preprocessing import StandardScaler

def convert_numpy(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj

class PhytoClust:
    def __init__(self, hmm_db, model_path, pathway_model_path, deep_model_path, cluster_max_distance=10000, min_genes_per_cluster=2):
        self.hmm_db = hmm_db
        self.cluster_model = joblib.load(model_path)
        self.pathway_model = joblib.load(pathway_model_path)
        self.deep_model = load_model(deep_model_path)
        self.scaler = StandardScaler()
        self.cluster_max_distance = cluster_max_distance
        self.min_genes_per_cluster = min_genes_per_cluster

        # Load real expression data
        self.expression_df = pd.read_csv("/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/metabolic/mgc_tools/phytoClust/arabidopsis_expression_mean_tpm.csv", index_col=0)

    def parse_fasta(self, fasta_file):
        return list(SeqIO.parse(fasta_file, "fasta"))

    def predict_genes(self, fasta_file):
        records = self.parse_fasta(fasta_file)
        predicted_genes = []
        for record in records:
            gene_id = record.id.split("|")[0] if "|" in record.id else record.id
            seq = str(record.seq)
            length = len(seq)
            gc = (seq.count("G") + seq.count("C")) / length if length else 0.0
            is_core = int(any(key in gene_id.upper() for key in ["CYP", "PKS", "TPS", "NRPS", "UGT", "ABC"]))

            # Use real expression if available
            expression = self.expression_df.loc[gene_id].values[0] if gene_id in self.expression_df.index else 0.0

            predicted_genes.append({
                "id": gene_id,
                "sequence": seq,
                "length": length,
                "gc_content": gc,
                "is_core": is_core,
                "expression": expression,
                "start": 0,
                "end": length
            })
        return predicted_genes

    def extract_kegg_proteins(self, kegg_path, output_fasta="kegg_proteins.fasta"):
        protein_records = []
        for species_dir in os.listdir(kegg_path):
            species_path = os.path.join(kegg_path, species_dir)
            if os.path.isdir(species_path):
                for filename in os.listdir(species_path):
                    if filename.endswith(".txt"):
                        with open(os.path.join(species_path, filename), "r") as f:
                            for line in f:
                                if not line.startswith("#"):
                                    parts = line.strip().split()
                                    if len(parts) > 1:
                                        gene_id, nucleotide_seq = parts[0], parts[1]
                                        protein_seq = str(Seq(nucleotide_seq).translate(to_stop=True))
                                        if len(protein_seq) > 50:
                                            protein_records.append(SeqRecord(Seq(protein_seq), id=gene_id, description=""))
        SeqIO.write(protein_records, output_fasta, "fasta")
        return output_fasta

    def run_hmmer_against_kegg(self, query_fasta, kegg_fasta="kegg_proteins.fasta", output_file="hmmer_kegg_results.txt"):
        cmd = f"phmmer --cpu 8 --tblout {output_file} {query_fasta} {kegg_fasta}"
        subprocess.run(cmd, shell=True, check=True)
        return output_file

    def parse_hmmer_results(self, hmmer_file, evalue_threshold=1e-5):
        metabolic_gene_scores = {}
        with open(hmmer_file, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    parts = line.split()
                    if len(parts) >= 6:
                        query_gene = parts[0]
                        e_value = float(parts[4])
                        bit_score = float(parts[5])
                        if e_value < evalue_threshold:
                            metabolic_gene_scores[query_gene] = bit_score
        return metabolic_gene_scores

    def detect_metabolic_gene_clusters(self, annotated_genes, metabolic_gene_scores):
        clusters = []
        metabolic_genes = [g for g in annotated_genes if g["id"] in metabolic_gene_scores]
        if not metabolic_genes:
            return clusters
        sorted_genes = sorted(metabolic_genes, key=lambda x: x["start"])
        current_cluster = []
        cluster_start = None
        for gene in sorted_genes:
            if not current_cluster:
                current_cluster.append(gene)
                cluster_start = gene["start"]
            else:
                if gene["start"] - cluster_start <= self.cluster_max_distance:
                    current_cluster.append(gene)
                else:
                    if len(current_cluster) >= self.min_genes_per_cluster:
                        clusters.append(current_cluster)
                    current_cluster = [gene]
                    cluster_start = gene["start"]
        if len(current_cluster) >= self.min_genes_per_cluster:
            clusters.append(current_cluster)
        return clusters

    def run_pipeline(self, fasta_file, kegg_path):
        kegg_fasta = self.extract_kegg_proteins(kegg_path)
        predicted_genes = self.predict_genes(fasta_file)
        SeqIO.write([SeqRecord(Seq(g['sequence']), id=g['id'], description="") for g in predicted_genes], "predicted_genes.fasta", "fasta")
        hmmer_results = self.run_hmmer_against_kegg("predicted_genes.fasta", kegg_fasta)
        metabolic_genes = self.parse_hmmer_results(hmmer_results)
        clusters = self.detect_metabolic_gene_clusters(predicted_genes, metabolic_genes)
        return clusters

if __name__ == "__main__":
    analyzer = PhytoClust(
        hmm_db="/groups/itay_mayrose/alongonda/datasets/hmmer_profile/Pfam-A.hmm",
        model_path="/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/metabolic/mgc_tools/phytoClust/models/phytoClust_optimized_model.pkl",
        pathway_model_path="/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/metabolic/mgc_tools/phytoClust/models/pathway_prediction_model.pkl",
        deep_model_path="/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/metabolic/mgc_tools/phytoClust/models/phytoClust_deep_model.h5"
    )

    clusters = analyzer.run_pipeline(
        fasta_file="/groups/itay_mayrose/alongonda/datasets/arabidopsis_thaliana/protein/Arabidopsis_thaliana.TAIR10.protein.fa",
        kegg_path="/groups/itay_mayrose/alongonda/datasets/KEGG"
    )

    with open("cluster_results.json", "w") as f:
        json.dump({"clusters": clusters}, f, indent=4, default=convert_numpy)

    logger.info("✅ Analysis complete! Results saved to cluster_results.json")
