import os
import pandas as pd
import matplotlib.pyplot as plt

def main():
    # Specify the root directory containing the CSV files
    root_directory = "/groups/itay_mayrose_nosnap/alongonda/Plant_MGC/sliding_window_outputs_with_statistics"
    file_pattern = "statistics"

    # Recursively find all CSV files matching the pattern
    csv_files = [
        os.path.join(dirpath, file)
        for dirpath, _, files in os.walk(root_directory)
        for file in files
        if file.startswith(file_pattern) and file.endswith(".csv")
    ]

    if not csv_files:
        print(f"No files matching pattern '{file_pattern}*.csv' found in {root_directory}")
        return

    # Consolidate data from all files
    rankings = []
    data = []
    for file in csv_files:
        df = pd.read_csv(file)
        if "enrichment_ratio" in df.columns:
            file_name = os.path.basename(file)
            median_enrichment = df["enrichment_ratio"].median()
            rankings.append({"file_name": file_name, "median_enrichment": median_enrichment})
            for enrichment_ratio in df["enrichment_ratio"]:
                data.append({"file_name": file_name, "enrichment_ratio": enrichment_ratio})

    # Create a DataFrame for ranking
    ranking_df = pd.DataFrame(rankings)
    ranking_df = ranking_df.sort_values(by="median_enrichment", ascending=False)

    # Save the ranking to a CSV file
    output_ranking_path = os.path.join(root_directory, "file_rankings.csv")
    ranking_df.to_csv(output_ranking_path, index=False)
    print(f"Ranking saved to {output_ranking_path}")

    # Create a consolidated DataFrame for plotting
    consolidated_df = pd.DataFrame(data)

    # Plot the enrichment_ratio vs file_name (ordered by ranking)
    ordered_files = ranking_df["file_name"]
    plt.figure(figsize=(12, 6))
    plt.boxplot(
        [
            consolidated_df[consolidated_df["file_name"] == name]["enrichment_ratio"].values
            for name in ordered_files
        ],
        labels=ordered_files,
        showfliers=False,
    )
    plt.xticks(rotation=90)
    plt.title("Enrichment Ratio by File Name (Ordered by Median)")
    plt.xlabel("File Name")
    plt.ylabel("Enrichment Ratio")
    plt.tight_layout()

    # Save the plot
    output_plot_path = output_ranking_path = os.path.join(root_directory, "enrichment_ratio_plot.png")
    plt.savefig(output_plot_path)
    print(f"Plot saved to {output_plot_path}")
    plt.show()

if __name__ == "__main__":
    main()
