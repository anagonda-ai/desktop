import os
import pandas as pd
import random
from typing import List
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# --- CONFIG --- #
WINDOW_SIZE = 10
WINDOWS_PER_FILE = 1000
MAX_WORKERS = 32
OUTPUT_PATH = "/groups/itay_mayrose/alongonda/Plant_MGC/kegg_metabolic_output/kegg_scanner_min_genes_based_metabolic/min_genes_3/random_control_windows.csv"
REQUIRED_COLUMNS = {"id", "start", "end", "chromosome", "metabolic_gene", "pathway"}

GENOME_DIRS = [
    "/groups/itay_mayrose/alongonda/datasets/full_genomes/ensembl/processed_annotations_test_no_chloroplast_with_sequences/metabolic_genes_with_letters",
    "/groups/itay_mayrose/alongonda/datasets/full_genomes/phytozome/processed_annotations_with_chromosomes_no_chloroplast_with_sequences/metabolic_genes_with_letters",
    "/groups/itay_mayrose/alongonda/datasets/full_genomes/plaza/processed_annotations_with_chromosomes_no_chloroplast_with_sequences/metabolic_genes_with_letters"
]

# --- DATA STRUCTURE --- #

@dataclass
class RandomWindow:
    genes: List[str]
    metabolic_genes: List[str]
    annotations: List[str]
    start: int
    end: int
    source_file: str

    def to_dict(self) -> dict:
        return {
            "genes": ",".join(self.genes),
            "metabolic_genes": ",".join(self.metabolic_genes),
            "metabolic_genes_annotations": ",".join(self.annotations),
            "start": self.start,
            "end": self.end,
            "source_file": self.source_file
        }

# --- SAMPLER CLASS --- #

class GenomeWindowSampler:
    def __init__(self, filepath: str, window_size: int, windows_per_file: int):
        self.filepath = filepath
        self.window_size = window_size
        self.windows_per_file = windows_per_file

    def is_valid(self, df: pd.DataFrame) -> bool:
        return REQUIRED_COLUMNS.issubset(df.columns) and not df.empty

    def sample_windows(self) -> List[RandomWindow]:
        try:
            df = pd.read_csv(self.filepath)
            if not self.is_valid(df):
                return []

            windows = []
            chromosomes = df["chromosome"].dropna().unique()
            if len(chromosomes) == 0:
                print(f"⚠️ No valid chromosomes in: {self.filepath}")
                return []
            windows_per_chrom = max(self.windows_per_file // len(chromosomes), 1)

            for chrom in chromosomes:
                chrom_df = df[df["chromosome"] == chrom].sort_values("start").reset_index(drop=True)
                if len(chrom_df) < self.window_size:
                    print(f"↪️ Skipping chromosome {chrom} in {self.filepath} (only {len(chrom_df)} genes)")
                    continue

                max_start = len(chrom_df) - self.window_size
                start_indices = random.sample(range(max_start), min(windows_per_chrom, max_start))

                for start_idx in start_indices:
                    window_df = chrom_df.iloc[start_idx:start_idx + self.window_size]
                    genes = window_df["id"].tolist()
                    # Use metabolic_gene column instead of annotation
                    metabolic_genes = window_df[window_df["metabolic_gene"].notna()]["id"].tolist()
                    annotations = window_df[window_df["metabolic_gene"].notna()]["pathway"].fillna("").astype(str).tolist()


                    window = RandomWindow(
                        genes=genes,
                        metabolic_genes=metabolic_genes,
                        annotations=annotations,
                        start=window_df["start"].min(),
                        end=window_df["end"].max(),
                        source_file=self.filepath
                    )
                    windows.append(window)

            return windows

        except Exception as e:
            print(f"❌ Error processing {self.filepath}: {e}")
            return []

# --- CONTROLLER --- #

class RandomWindowGenerator:
    def __init__(self, genome_dirs: List[str]):
        self.genome_files = self.collect_genome_files(genome_dirs)

    def collect_genome_files(self, dirs: List[str]) -> List[str]:
        all_files = []
        for d in dirs:
            all_files += [os.path.join(d, f) for f in os.listdir(d) if f.endswith(".csv")]
        return all_files

    def generate_all(self) -> pd.DataFrame:
        print(f"🧬 Found {len(self.genome_files)} genome files.")
        windows = []

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {
                executor.submit(self._process_file, f): f for f in self.genome_files
            }
            for future in tqdm(as_completed(futures), total=len(futures), desc="Sampling random windows", unit="file"):
                result = future.result()
                if result:
                    windows.extend(result)

        df = pd.DataFrame([w.to_dict() for w in windows])
        df.drop_duplicates(inplace=True)
        df.sort_values(by=["source_file", "start"], inplace=True)
        return df

    def _process_file(self, filepath: str) -> List[RandomWindow]:
        sampler = GenomeWindowSampler(filepath, WINDOW_SIZE, WINDOWS_PER_FILE)
        return sampler.sample_windows()

# --- MAIN --- #

def main():
    generator = RandomWindowGenerator(GENOME_DIRS)
    df = generator.generate_all()

    if not df.empty:
        df.to_csv(OUTPUT_PATH, index=False)
        print(f"✅ Saved {len(df)} random windows → {OUTPUT_PATH}")
    else:
        print("⚠️ No random windows were generated.")

if __name__ == "__main__":
    main()
