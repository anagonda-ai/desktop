import csv


def filter_blast_results(input_file, output_file):
    """Filters BLAST results to keep only perfect hits (100% identity and coverage).

    Args:
        input_file (str): The path to the input BLAST results file.
        output_file (str): The path to the output file for filtered results.
    """

    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        reader = csv.reader(f_in, delimiter='\t')
        writer = csv.writer(f_out, delimiter='\t')

        header = next(reader)
        writer.writerow(header)


        for row in reader:
            query_acc, subject_acc, identity, alignment_length, mismatches, gaps, qstart, qend, sstart, send, evalue, bitscore = row

            # Check for perfect hit (100% identity and coverage)
            if float(identity) == 100 and int(alignment_length) == int(qend) - int(qstart) + 1:
                writer.writerow(row)

# Replace with your input and output file paths
input_file = "/groups/itay_mayrose/alongonda/datasets/plantcyc/pmn_on_mibig.fasta"
output_file = "/groups/itay_mayrose/alongonda/datasets/plantcyc/pmn_on_mibig_perfect_hits.fasta"

filter_blast_results(input_file, output_file)