import os, glob, csv
import concurrent.futures

BLAST_ROOT = "/groups/itay_mayrose/alongonda/datasets/KEGG_annotations_modules_metabolic/selected_genes/aans_selected_1000/best_hits/"


def process_result_file(args):
    qbase, fpath = args
    local_best = {}
    try:
        with open(fpath, "r") as fh:
            for ln in fh:
                cols = ln.rstrip().split("\t")
                if len(cols) < 12:
                    continue

                qid    = cols[0]
                sid    = cols[1]
                pident = cols[2]
                length = cols[3]
                mismatch = cols[4]
                gapopen  = cols[5]
                qstart   = cols[6]
                qend     = cols[7]
                sstart   = cols[8]
                send     = cols[9]
                evalue   = float(cols[10])
                bitscore = float(cols[11])

                key = (qbase, qid, fpath)
                rec = dict(
                    query_dir=qbase, query_gene=qid, subject_gene=sid,
                    pident=pident, length=length, mismatch=mismatch, gapopen=gapopen,
                    qstart=qstart, qend=qend, sstart=sstart, send=send,
                    evalue=evalue, bitscore=bitscore, src=fpath
                )

                if key not in local_best or \
                   bitscore > local_best[key]['bitscore'] or \
                   (bitscore == local_best[key]['bitscore'] and evalue < local_best[key]['evalue']):
                    local_best[key] = rec
    except Exception as e:
        print("WARN reading", fpath, e)
    return local_best

# Gather all (qbase, fpath) pairs for BLAST result files
file_args = []
for qdir in sorted(glob.glob(os.path.join(BLAST_ROOT, "*"))):
    if not os.path.isdir(qdir):
        continue
    qbase = os.path.basename(qdir)
    for fpath in sorted(glob.glob(os.path.join(qdir, "*_results.txt"))):
        file_args.append((qbase, fpath))

# Use ThreadPoolExecutor for concurrency (since this is I/O bound)
best_per_query = {}
with concurrent.futures.ThreadPoolExecutor() as executor:
    futures = [executor.submit(process_result_file, arg) for arg in file_args]
    for fut in concurrent.futures.as_completed(futures):
        local_best = fut.result()
        for key, rec in local_best.items():
            if key not in best_per_query or \
               rec['bitscore'] > best_per_query[key]['bitscore'] or \
               (rec['bitscore'] == best_per_query[key]['bitscore'] and rec['evalue'] < best_per_query[key]['evalue']):
                best_per_query[key] = rec

# Create a mapping: organism -> list of records
from collections import defaultdict
import re

organism_records = defaultdict(list)

for key, rec in best_per_query.items():
    source_file = os.path.basename(rec["src"])
    organism = "_".join(source_file.split("_")[2:])
    organism_records[organism].append(rec)

# Write separate output file for each organism
for organism, records in organism_records.items():
    out_path = os.path.join(BLAST_ROOT, f"best_hits_{organism}.csv")
    with open(out_path, "w", newline='', encoding="utf-8") as outcsv:
        writer = csv.writer(outcsv)
        writer.writerow([
            "query_dir","query_gene","subject_gene","pident","length",
            "mismatch","gapopen","qstart","qend","sstart","send",
            "evalue","bitscore","source_file"
        ])
        for rec in records:
            writer.writerow([
                rec['query_dir'], rec['query_gene'], rec['subject_gene'],
                rec['pident'], rec['length'], rec['mismatch'], rec['gapopen'],
                rec['qstart'], rec['qend'], rec['sstart'], rec['send'],
                rec['evalue'], rec['bitscore'], rec['src']
            ])
    print(f"Wrote per-organism file: {out_path}")
