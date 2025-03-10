import re
import os

def convert_to_fasta(input_file, output_file, line_length=60):
    with open(input_file, 'r', errors='ignore') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            match = re.match(r'(\S+):\s+(.+)', line)
            if match:
                identifier = match.group(1)
                sequence = match.group(2)
                outfile.write(f'>{identifier}\n')
                # Split the sequence into lines of specified length
                for i in range(0, len(sequence), line_length):
                    outfile.write(sequence[i:i+line_length] + '\n')

def get_directory_path(file_path):
    return os.path.dirname(file_path)

def main():
    input_file = input("Enter the path to the input file: ")
    output_file = os.path.join(get_directory_path(input_file), 'output_file.fasta')
    line_length = 60  # Adjust this value as needed (e.g., 80 or None for no line breaks)
    convert_to_fasta(input_file, output_file, line_length)
    print(f"Output written to: {output_file}")
    
if __name__ == '__main__':
    main()