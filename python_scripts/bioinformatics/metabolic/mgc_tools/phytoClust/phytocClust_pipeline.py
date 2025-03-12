import os
import json
import subprocess
import numpy as np
import joblib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tensorflow.keras.models import load_model
from Bio import SeqIO
from sklearn.preprocessing import StandardScaler

# Convert NumPy data types to regular Python types
def convert_numpy(obj):
    if isinstance(obj, np.integer):
        return int(obj)  # Convert int64 → int
    elif isinstance(obj, np.floating):
        return float(obj)  # Convert float64 → float
    elif isinstance(obj, np.ndarray):
        return obj.tolist()  # Convert NumPy array → Python list
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

    def parse_fasta(self, fasta_file):
        return list(SeqIO.parse(fasta_file, "fasta"))

    def predict_genes(self, fasta_file):
        records = self.parse_fasta(fasta_file)
        predicted_genes = []
        for record in records:
            for strand, start in [(+1, 0), (-1, 1)]:
                for i in range(start, len(record.seq), 3):
                    codon = record.seq[i:i+3]
                    if codon == "ATG":
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
    
    def extract_kegg_proteins(kegg_path, output_fasta="kegg_proteins.fasta"):
        """Extract and translate KEGG gene sequences to protein sequences."""
        protein_records = []

        for species_dir in os.listdir(kegg_path):
            species_path = os.path.join(kegg_path, species_dir)
            if os.path.isdir(species_path):
                for filename in os.listdir(species_path):
                    if filename.endswith(".txt"):
                        file_path = os.path.join(species_path, filename)
                        with open(file_path, "r") as f:
                            for line in f:
                                if not line.startswith("#"):
                                    parts = line.strip().split()
                                    if len(parts) > 1:
                                        gene_id = parts[0]
                                        nucleotide_seq = parts[1]
                                        protein_seq = str(Seq(nucleotide_seq).translate(to_stop=True))
                                        
                                        if len(protein_seq) > 50:  # Ignore short sequences
                                            protein_records.append(
                                                SeqRecord(Seq(protein_seq), id=gene_id, description="")
                                            )

        SeqIO.write(protein_records, output_fasta, "fasta")
        print(f"✅ Extracted {len(protein_records)} protein sequences from KEGG.")

        return output_fasta
    
    def run_hmmer_against_kegg(query_fasta, kegg_fasta="kegg_proteins.fasta", output_file="hmmer_kegg_results.txt"):
        """Run HMMER to match predicted genes against KEGG metabolic proteins."""
        cmd = f"hmmsearch --cpu 8 --tblout {output_file} {query_fasta} {kegg_fasta}"
        subprocess.run(cmd, shell=True)
        return output_file

    def annotate_genes(self, genes, fasta_file):
        """Uses HMMER to scan genes against the Pfam-A.hmm database."""
        output_file = "hmmscan_results.txt"
        cmd = f"hmmsearch --cpu 8 --tblout {output_file} {self.hmm_db} {fasta_file}"
        subprocess.run(cmd, shell=True)

        annotations = {}
        with open(output_file, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    parts = line.split()
                    annotations[parts[3]] = parts[0]  # Query gene ID → Matched Pfam domain
        for gene in genes:
            gene["annotation"] = annotations.get(gene["id"], "Unknown")
        return genes

    def predict_pathways(self, genes):
        """Batch prediction for pathway classification using deep learning."""
        # Prepare batch feature matrix
        feature_matrix = np.array([[len(g['sequence']), 1, 0.42, 8.5] for g in genes])

        # Ensure correct shape (batch_size, num_features)
        feature_matrix = feature_matrix.reshape(len(genes), 4)

        # Run batch prediction once
        pathway_predictions = np.argmax(self.deep_model.predict(feature_matrix), axis=1)

        # Assign predicted pathways to genes
        for gene, predicted_pathway in zip(genes, pathway_predictions):
            gene['pathway'] = predicted_pathway

        return genes
    
    def parse_hmmer_results(hmmer_file, evalue_threshold=1e-5):
        """Parse HMMER results and return genes with strong KEGG similarity."""
        metabolic_genes = set()

        with open(hmmer_file, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    parts = line.split()
                    query_gene = parts[0]  # Our predicted gene
                    e_value = float(parts[6])  # E-value

                    if e_value < evalue_threshold:  # Strong match
                        metabolic_genes.add(query_gene)

        print(f"✅ Found {len(metabolic_genes)} metabolic genes based on KEGG similarity.")
        return metabolic_genes


    def detect_metabolic_gene_clusters(self, annotated_genes, metabolic_gene_set):
        """Detects clusters of metabolic genes within a 10kbp window."""
        clusters = []

        # Keep only genes that matched KEGG via HMMER
        metabolic_genes = [g for g in annotated_genes if g["id"] in metabolic_gene_set]

        if not metabolic_genes:
            print("⚠️ No metabolic genes found based on sequence similarity.")
            return clusters  

        # Sort genes by genomic position
        sorted_genes = sorted(metabolic_genes, key=lambda x: x['start'])

        current_cluster = []
        cluster_start = None

        for gene in sorted_genes:
            if not current_cluster:
                current_cluster.append(gene)
                cluster_start = gene['start']
            else:
                last_gene = current_cluster[-1]
                
                # Check if within 10kbp
                if gene['start'] - last_gene['end'] <= 10000:
                    current_cluster.append(gene)
                else:
                    # Save valid clusters
                    if len(current_cluster) >= self.min_genes_per_cluster:
                        clusters.append(current_cluster)
                    # Start a new cluster
                    current_cluster = [gene]
                    cluster_start = gene['start']

        # Add last cluster if it meets the criteria
        if len(current_cluster) >= self.min_genes_per_cluster:
            clusters.append(current_cluster)

        return clusters



    def run_pipeline(self, fasta_file, kegg_path):
        """Complete pipeline: predict genes, match them to KEGG, and cluster metabolic genes."""
        kegg_fasta = self.extract_kegg_proteins(kegg_path)  # Convert KEGG sequences to protein
        predicted_genes = self.predict_genes(fasta_file)  # Predict genes from FASTA
        predicted_fasta = "predicted_genes.fasta"
        
        # Write predicted genes to FASTA for HMMER
        SeqIO.write([SeqRecord(Seq(g['sequence']), id=g['id'], description="") for g in predicted_genes], 
                    predicted_fasta, "fasta")

        hmmer_results = self.run_hmmer_against_kegg(predicted_fasta, kegg_fasta)  # Match genes to KEGG
        metabolic_genes = self.parse_hmmer_results(hmmer_results)  # Extract metabolic genes
        
        annotated_genes = self.annotate_genes(predicted_genes, fasta_file, metabolic_genes)
        clusters = self.detect_metabolic_gene_clusters(annotated_genes, metabolic_genes)

        return clusters

if __name__ == "__main__":
    hmm_db = "/groups/itay_mayrose/alongonda/datasets/hmmer_profile/Pfam-A.hmm"
    analyzer = PhytoClust(hmm_db, "/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/metabolic/mgc_tools/phytoClust/models/phytoClust_optimized_model.pkl", "/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/metabolic/mgc_tools/phytoClust/models/pathway_prediction_model.pkl", "/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/metabolic/mgc_tools/phytoClust/models/phytoClust_deep_model.h5")

    genome_file = "/groups/itay_mayrose/alongonda/datasets/arabidopsis_thaliana/protein/Arabidopsis_thaliana.TAIR10.protein.fa"
    clusters = analyzer.run_pipeline(genome_file)

    with open("cluster_results.json", "w") as f:
        json.dump({"clusters": clusters}, f, indent=4, default=convert_numpy)

    print("✅ Analysis complete! Results saved to cluster_results.json")
