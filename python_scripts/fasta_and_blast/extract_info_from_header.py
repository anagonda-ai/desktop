import concurrent.futures
import csv
import os

def extract_info_from_header(header):
    # Regular expression to extract gene and transcript name, start, and end positions
    parts = header.split(' ')
    gene_name = parts[3].split(':')[1]
    transcript_name = parts[4].split(':')[1]
    start_position = int(parts[2].split(':')[3])
    end_position = int(parts[2].split(':')[4])
    return gene_name, transcript_name, start_position, end_position

def process_file(file_path):
    data = []
    with open(file_path, 'r') as f:
        # Variables to store current entry information
        current_header = None
        current_sequence = []

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # If we already have a previous entry, save it
                if current_header and current_sequence:
                    gene_name, transcript_name, start_position, end_position = extract_info_from_header(current_header)
                    protein_sequence = ''.join(current_sequence)
                    
                    if gene_name and transcript_name and start_position and end_position:
                        data.append((gene_name, transcript_name, start_position, end_position, protein_sequence))
                
                # Reset for new entry
                current_header = line
                current_sequence = []
            else:
                # Collect sequence lines
                current_sequence.append(line)
        
        # Process the last entry
        if current_header and current_sequence:
            gene_name, transcript_name, start_position, end_position = extract_info_from_header(current_header)
            protein_sequence = ''.join(current_sequence)
            
            if gene_name and transcript_name and start_position and end_position:
                data.append((gene_name, transcript_name, start_position, end_position, protein_sequence))
    
    # Sort the data by start and end values
    data.sort(key=lambda x: (x[2], x[3]))
    
    # Write the data to a CSV file
    output_file = os.path.join(os.path.dirname(file_path), 'extracted_data.csv')
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['gene_name', 'transcript', 'start', 'end', 'protein_sequence'])
        writer.writerows(data)
    
    print(f"Data saved to {output_file}")

def process_fa_files(root_dir):
    file_paths = []
    for subdir, _, files in os.walk(root_dir):
        for file in files:
            if file.endswith('.fa'):
                file_path = os.path.join(subdir, file)
                file_paths.append(file_path)
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = {executor.submit(process_file, file_path): file_path for file_path in file_paths}
        for future in concurrent.futures.as_completed(futures):
            file_path = futures[future]
            try:
                future.result()
                print(f"Processed file: {file_path}")
            except Exception as exc:
                print(f"Error processing file {file_path}: {exc}")

def main():
    root_dir = "/groups/itay_mayrose_nosnap/alongonda/full_genomes/ensembl/organisms"
    process_fa_files(root_dir)

if __name__ == "__main__":
    main()