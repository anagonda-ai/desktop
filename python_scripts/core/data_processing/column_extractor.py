import csv
import re
import pandas as pd
import os
import sys

def remove_parentheses_content(text):
    # Use regular expression to remove content between parentheses
    return re.sub(r'\(.*?\)', '', text)

def process_csv_file(csv_file_path, output_csv_file_path):
    data = []
    # Increase the CSV field size limit
    max_int = sys.maxsize
    while True:
        try:
            csv.field_size_limit(max_int)
            break
        except OverflowError:
            max_int = int(max_int / 10)

    with open(csv_file_path, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        next(csv_reader)  # Skip the header row
        for row in csv_reader:
            if len(row) >= 3:
                row_str = '\t'.join(row)
                row_str = remove_parentheses_content(row_str).split('\t')
                id = row_str[0]
                start = row_str[3]
                end = row_str[4]
                sequence = row_str[5]
                data.append([id, start, end, sequence])

    # Create a new DataFrame
    df = pd.DataFrame(data, columns=['id', 'start', 'end', 'sequence'])
    df.to_csv(output_csv_file_path, index=False)
    print(f"Data saved to {output_csv_file_path}")

def main():
    directory = "/groups/itay_mayrose/alongonda/datasets/full_genomes/ensembl/organisms"
    output_directory = "/groups/itay_mayrose/alongonda/datasets/full_genomes/ensembl/processed_annotations_test"
    os.makedirs(output_directory, exist_ok=True)

    for root, dirs, files in os.walk(directory):
        for filename in files:
            if filename.endswith(".csv"):
                csv_file_path = os.path.join(root, filename)
                output_file_name = os.path.basename(os.path.dirname(csv_file_path))
                output_csv_file_path = os.path.join(output_directory, output_file_name + ".csv")
                process_csv_file(csv_file_path, output_csv_file_path)

if __name__ == "__main__":
    main()