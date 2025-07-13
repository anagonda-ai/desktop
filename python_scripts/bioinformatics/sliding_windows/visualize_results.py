import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Path to the results CSV file
csv_file = '/groups/itay_mayrose/alongonda/Plant_MGC/kegg_verified_scanner_min_genes_3_overlap_merge/kegg_scanner_min_genes_based_metabolic/min_genes_3/potential_groups_w10_filtered.csv'

# Load the CSV file
df = pd.read_csv(csv_file)

def clean_filename(raw_path):
    name = os.path.basename(raw_path).lower()
    for suffix in [
        '_updated_annotated.csv', '_annotated.csv', '.filtered.csv',
        '.filtered', '.csv'
    ]:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    name = name.strip('._')
    name = name.replace('..', '.')
    name = name.replace('_', ' ')
    name = name.title()
    return name

if 'source_file' in df.columns:
    counts = df['source_file'].value_counts()
    normalized_counts = {}
    normalized_data = []

    for source_file, count in counts.items():
        source_path = source_file if os.path.isabs(source_file) else os.path.join(os.path.dirname(csv_file), source_file)
        try:
            file_df = pd.read_csv(source_path)
            file_rows = len(file_df)
            normalized = count / file_rows if file_rows > 0 else 0
        except Exception as e:
            print(f"Could not read {source_path}: {e}")
            file_rows = None
            normalized = None

        cleaned_name = clean_filename(source_file)
        normalized_counts[cleaned_name] = normalized
        print(f"Normalized count for {cleaned_name}: {normalized}")

        if file_rows is not None and normalized is not None:
            normalized_data.append((cleaned_name, file_rows, normalized))

    # Build series and DataFrame
    normalized_counts = pd.Series(normalized_counts)
    df_corr = pd.DataFrame(normalized_data, columns=["source", "file_rows", "normalized"])

    # Plot normalized bar chart
    plt.figure(figsize=(10, 20))
    normalized_counts_sorted = normalized_counts.sort_values(ascending=True)
    normalized_counts_sorted.plot(kind='barh', color='seagreen', edgecolor='black')
    plt.title('Normalized Results per Source File (Sorted)')
    plt.xlabel('Normalized Number of Results')
    plt.ylabel('Source File')
    plt.tight_layout()
    bar_output_path = os.path.join(os.path.dirname(csv_file), 'results_per_source_file.png')
    plt.savefig(bar_output_path)
    plt.show()
    print(f"✅ Bar figure saved to: {bar_output_path}")

    # Compute and print correlation
    corr_value, p_value = pearsonr(df_corr["file_rows"], df_corr["normalized"])
    print(f"\n📊 Pearson correlation between file size (rows) and normalized count: r = {corr_value:.3f}, p = {p_value:.2e}")

    # Optional scatter plot
    plt.figure(figsize=(6, 5))
    plt.scatter(df_corr["file_rows"], df_corr["normalized"], alpha=0.6)
    plt.title(f"Correlation: File Rows vs. Normalized Count\nr={corr_value:.3f}, p={p_value:.1e}")
    plt.xlabel("Number of Rows in Source File")
    plt.ylabel("Normalized Result Count")
    plt.grid(True)
    plt.tight_layout()
    scatter_output_path = os.path.join(os.path.dirname(csv_file), 'correlation_scatter.png')
    plt.savefig(scatter_output_path)
    plt.show()
    print(f"✅ Scatter figure saved to: {scatter_output_path}")

else:
    print("❌ No 'source_file' column found in the CSV.")
