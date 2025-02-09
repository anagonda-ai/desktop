import dask.array as da
import dask.dataframe as dd
import pandas as pd
import numpy as np
import hdbscan
import multiprocessing
from sklearn.preprocessing import MinMaxScaler

# ----------------- Load Data Efficiently ----------------- #
def load_data(csv_file):
    """Load CSV data using Dask for memory efficiency."""
    print("Loading data...")
    ddf = dd.read_csv(csv_file, usecols=["Similarity"], dtype={'Similarity': 'float32'})
    return ddf

# ----------------- Normalize Data ----------------- #
def normalize_data(df):
    """Normalize the similarity scores using MinMaxScaler."""
    print("Normalizing data...")
    scaler = MinMaxScaler()

    # Convert the 'Similarity' column to a Dask Array
    similarity_array = df["Similarity"].to_dask_array(lengths=True)

    # Fit the scaler on the entire data (compute necessary to fit on all data)
    similarity_array_computed = similarity_array.compute()
    scaler.fit(similarity_array_computed.reshape(-1, 1))

    # Transform the data using the fitted scaler
    normalized_similarity = scaler.transform(similarity_array_computed.reshape(-1, 1)).squeeze()

    # Create a Dask Array from the normalized data
    normalized_similarity_dask = da.from_array(normalized_similarity, chunks=similarity_array.chunks)

    # Ensure the Dask Array has the same number of partitions as the DataFrame
    normalized_similarity_dask = normalized_similarity_dask.rechunk(df.npartitions)

    # Assign the normalized similarity scores back to the DataFrame
    df["normalized_similarity"] = normalized_similarity_dask

    return df

# ----------------- Perform Clustering Using HDBSCAN ----------------- #
def perform_clustering(df, min_cluster_size=10, min_samples=5):
    """Perform HDBSCAN clustering on the entire dataset."""
    print("Performing HDBSCAN clustering...")

    # Convert Dask dataframe to Pandas (HDBSCAN does not support Dask)
    df = df.compute()  

    # Prepare feature array for HDBSCAN
    X = df[["normalized_similarity"]].values

    # Run HDBSCAN clustering
    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples, metric='euclidean', core_dist_n_jobs=multiprocessing.cpu_count())
    df["cluster"] = clusterer.fit_predict(X)

    return df

# ----------------- Save Results Efficiently ----------------- #
def save_results(df, output_csv):
    """Save clustered results to a CSV file."""
    print("Saving results...")
    df.to_csv(output_csv, index=False)
    print(f"Clustered data saved to {output_csv}")

# ----------------- Main Function ----------------- #
def main(csv_file, output_csv, min_cluster_size=10, min_samples=5):
    """Load, process, cluster, and save the entire dataset."""
    
    # Load Data
    ddf = load_data(csv_file)
    
    # Normalize Data
    df_normalized = normalize_data(ddf)
    
    # Perform Clustering
    df_clustered = perform_clustering(df_normalized, min_cluster_size, min_samples)
    
    # Save Results
    save_results(df_clustered, output_csv)

if __name__ == "__main__":
    csv_file = '/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/protein_global_alignment/graph_output.csv'
    output_csv = '/groups/itay_mayrose/alongonda/desktop/python_scripts/bioinformatics/protein_global_alignment/ml_clustering_output/clustered_output.csv'

    main(csv_file, output_csv)
