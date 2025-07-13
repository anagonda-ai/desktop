import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm


def fasta_to_dataframe(fasta_path):
    """
    Load a FASTA file into a pandas DataFrame.
    Supports two header formats:
    1. mgc_candidate | gene_id | kegg_id | source_file
    2. gene_id | mgc_id | start | end
    """
    records = []
    with open(fasta_path, 'r') as f:
        header = None
        seq_lines = []

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    records.append(parse_record(header, seq_lines))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)

        if header:
            records.append(parse_record(header, seq_lines))

    df = pd.DataFrame(records)
    return df


def parse_record(header, seq_lines):
    """
    Parse header line and sequence into structured fields.
    Supports two formats:
    - mgc_candidate | gene_id | kegg_id | source_file
    - gene_id | mgc_id | start | end
    """
    parts = [part.strip() for part in header.split('|')]
    sequence = ''.join(seq_lines)

    if len(parts) != 4:
        raise ValueError(f"Invalid header format: {header}")

    # Try to detect if start and end are integers
    try:
        start = int(parts[2])
        end = int(parts[3])
        # Format is: gene_id | mgc_id | start | end
        return {
            'gene_id': parts[0],
            'mgc_id': parts[1],
            'start': start,
            'end': end,
            'sequence': sequence
        }
    except ValueError:
        # Not integers → format is: mgc_candidate | gene_id | kegg_id | source_file
        return {
            'mgc_candidate': parts[0],
            'gene_id': parts[1],
            'kegg_id': parts[2],
            'source_file': parts[3],
            'sequence': sequence
        }


def process_fasta_file(fasta_file, output_dir):
    """
    Convert a single FASTA file to CSV in the output directory.
    """
    try:
        df = fasta_to_dataframe(fasta_file)
        filename = os.path.splitext(os.path.basename(fasta_file))[0] + '.csv'
        output_path = os.path.join(output_dir, filename)
        df.to_csv(output_path, index=False)
        return output_path
    except Exception as e:
        print(f"[✘] Failed: {os.path.basename(fasta_file)} - {e}")
        return None


def process_fasta_list(list_txt_path, output_dir, max_workers=8):
    """
    Process a list of FASTA file paths from a text file and convert them to CSV.
    """
    if not os.path.isfile(list_txt_path):
        raise FileNotFoundError(f"List file not found: {list_txt_path}")

    with open(list_txt_path, 'r') as f:
        fasta_files = [line.strip() for line in f if line.strip()]

    print(f"[INFO] Found {len(fasta_files)} FASTA files to process.")

    os.makedirs(output_dir, exist_ok=True)

    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_fasta_file, fasta, output_dir): fasta
            for fasta in fasta_files
        }

        with tqdm(total=len(futures), desc="Processing FASTA files", ncols=100) as pbar:
            for future in as_completed(futures):
                fasta = futures[future]
                try:
                    result = future.result()
                    if result:
                        results.append(result)
                except Exception as e:
                    print(f"[✘] Error processing {fasta}: {e}")
                finally:
                    pbar.update(1)

    print(f"[DONE] Processed {len(results)} files successfully.")


if __name__ == "__main__":
    fasta_list_path = '/groups/itay_mayrose/alongonda/Plant_MGC/kegg_metabolic_output_g3_slurm_no_chloroplast/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/merged_list.txt'

    output_folder = '/groups/itay_mayrose/alongonda/Plant_MGC/kegg_metabolic_output_g3_slurm_no_chloroplast/kegg_scanner_min_genes_based_metabolic/min_genes_3/mgc_candidates_fasta_files_without_e2p2_filtered_test/mgc_candidates_csv_files'

    process_fasta_list(fasta_list_path, output_folder, max_workers=30)
