def merge_fasta_files(files, output_file):
    with open(output_file, 'w') as outfile:
        for file in files:
            with open(file, 'r') as infile:
                outfile.write(infile.read())

def main():
    files = []
    while True:
        file_path = input("Enter the path to a FASTA file (or -1 to finish): ")
        if file_path == '-1':
            break
        files.append(file_path)
    
    output_file = input("Enter the path to the output FASTA file: ")
    
    merge_fasta_files(files, output_file)
    print(f"Merged FASTA file written to: {output_file}")

if __name__ == '__main__':
    main()