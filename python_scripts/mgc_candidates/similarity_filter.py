import logging
import os
import argparse
import itertools
import subprocess
from Bio import SeqIO, SearchIO
from pathlib import Path
from collections import defaultdict
import concurrent.futures

def locate_filters(filter_by, main_path, output):
    temp_fb = []
    for fb in filter_by:
        if fb.is_dir():
            temp = main_path / output / f'{fb.name}.fasta'
            with open(temp, 'w') as out_file:
                for file in fb.iterdir():
                    records = list(SeqIO.parse(file, 'fasta'))
                    for record in records:
                        record.id = f'{record.id}_{file.stem}'
                        record.description = ''
                    SeqIO.write(records, out_file, 'fasta')
                print(f'Combined {fb.stem} filter-by files into {temp.name}')
            temp_fb.append(temp)
        elif fb.is_file():
            print(f'Found {fb.name} as a filter-by file')
            temp_fb.append(fb)
    return temp_fb

def create_blast_db(pt_genes_file, db_name):
    result = subprocess.run(["makeblastdb", "-in", pt_genes_file, "-dbtype", "prot", "-out", db_name], check=True)
    if result.returncode != 0:
        print("makeblastdb error:", result.stderr)

def run_blast(targets_file, db_name, output_file):
    try:
        result = subprocess.run([
            "blastp",
            "-query", targets_file,
            "-db", db_name,
            "-out", output_file,
            "-outfmt", "5"  # Explicitly specify XML output format
        ], capture_output=True, text=True, check=True)
        logging.info(f'Completed BLAST search: {output_file}')
        
        # Verify the output file exists and has content
        if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            raise subprocess.SubprocessError(f"BLAST output file {output_file} is empty or not created")
            
    except subprocess.CalledProcessError as e:
        logging.error(f'blastp error: {e.stderr}')
        raise

def parse_blast_results(blast_output, id_map, evalue_threshold=0.001):
    hits = set()
    try:
        # Use parse instead of read for multiple queries
        for blast_qresult in SearchIO.parse(blast_output, 'blast-xml'):
            try:
                query_id = blast_qresult.id
                # If the query_id is directly in id_map, use it
                if query_id in id_map:
                    current_id = query_id
                else:
                    # Try to parse index if direct lookup fails
                    try:
                        query_index = int(query_id.split('_')[1]) - 1
                        current_id = id_map[query_index]
                    except (IndexError, ValueError):
                        print(f'Could not parse query ID: {query_id}')
                        continue

                for hit in blast_qresult:
                    for hsp in hit.hsps:
                        if hsp.evalue < evalue_threshold:
                            hits.add(current_id)
                            break
            except Exception as e:
                print(f'Error processing hit in query {blast_qresult.id}: {str(e)}')
                continue
    except Exception as e:
        print(f'Error parsing BLAST results from {blast_output}: {str(e)}')
        # Try text format parsing as fallback
        try:
            with open(blast_output, 'r') as f:
                current_query = None
                for line in f:
                    if line.startswith('Query='):
                        # Extract the query ID
                        query_id = line.split('=')[1].strip().split()[0]
                        if query_id in id_map:
                            current_query = query_id
                        else:
                            try:
                                query_index = int(query_id.split('_')[1]) - 1
                                current_query = id_map[query_index]
                            except (IndexError, ValueError):
                                current_query = None
                                print(f'Could not parse query ID: {query_id}')
                    elif line.strip() and not line.startswith(('#', 'Database:', 'Length=', 'BLASTP', ' ')):
                        try:
                            fields = line.split()
                            if len(fields) >= 3 and float(fields[-1]) < evalue_threshold:
                                if current_query:
                                    hits.add(current_query)
                        except (ValueError, IndexError):
                            continue
        except Exception as text_error:
            print(f'Text parsing also failed: {str(text_error)}')
            raise
    return hits

def split_by_combinations(input_dict):
    all_keys = list(input_dict.keys())
    result = defaultdict(list)
    
    # iterate over all possible combinations of the keys:
    # reverse order is important!!! (because we remove elements from the lists)
    for r in range(len(all_keys), 0, -1):
        for combo in itertools.combinations(all_keys, r):

            intersected_values = set(input_dict[combo[0]])
            for key in combo[1:]:
                intersected_values &= set(input_dict[key])

            if intersected_values:
                result[' & '.join(combo)].extend(intersected_values)
                # remove the intersected values from the original lists to ensure uniqueness
                for key in combo:
                    input_dict[key] = list(set(input_dict[key]) - intersected_values)
    
    return dict(result)

def process_target(target, output, main_path, filter_by, silent):
    target_name = os.path.basename(target)
    os.mkdir(f"{output}/{target_name}")
    blast_db = main_path/output/target_name/'blast_db'  
    blast_output = main_path/output/target_name/f'blast_results_{target_name}.xml'
    output_file = main_path/output/target_name/f'blast_filtered_targets_{target_name}.fasta'
    id_map = [record.id for record in SeqIO.parse(target, 'fasta')]
    hits = {}
    for fb in filter_by:
        create_blast_db(fb, blast_db)
        run_blast(target, blast_db, blast_output)
        hits[fb.stem] = parse_blast_results(blast_output, id_map)

    all_hits = set.union(*hits.values())
    hits = split_by_combinations(hits)

    print(f'Filtered {len(all_hits)} out of {len(id_map)} sequences:')
    if silent:
        for k, v in hits.items():
            print(f'\tFrom {k}: {len(v)} hits')
    else:
        for k, v in hits.items():
            print(f"\tFrom {k}:\n\t\t{'    '.join(v)}")
    
    filtered_records = [record for record in SeqIO.parse(target, 'fasta') if record.id not in all_hits]
    SeqIO.write(filtered_records, output_file, 'fasta')
    print(f'Filtered sequences written to {output_file.name}')
    
def filter_by_similarity(targets_dir, filter_by, relative, output, silent=True):
    
    os.mkdir(output)

    main_path = Path(__file__).parent
    if not isinstance(filter_by, list):
        filter_by = [filter_by]

    if relative:
        targets_dir = main_path / targets_dir
        filter_by = [main_path / f for f in filter_by]
    else:
        targets_dir = Path(targets_dir)
        filter_by = [Path(f) for f in filter_by]

    filter_by = locate_filters(filter_by, main_path, output)

      
    target_fastas = list(targets_dir.rglob('*.fasta'))
    target_fastas = [f for f in target_fastas if f.name.startswith("start_end_")]
    
    print(f'Found {len(target_fastas)} target files in {targets_dir}')
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_target, target, output, main_path, filter_by, silent) for target in target_fastas]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as exc:
                print(f'Generated an exception: {exc}')

def main(args):

    filter_by_similarity(args.targets_dir, args.filter_by, args.relative, args.output, silent=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter by similarity')
    parser.add_argument('-t', '--targets_dir', type=str, default='targets.fa', help='targets fasta dir to be filtered')
    parser.add_argument('-f', '--filter_by', type=str, nargs='+', default='pt_genes', \
                        help='list of: fasta-containing folder or combined file with the filter-by sequences')
    parser.add_argument('-o', '--output', type=str, default='', help='outputs prefix')
    parser.add_argument('--relative', action='store_true', default=True, \
                        help='use flag if all required inputs are in the script directory')    
    main(parser.parse_args())