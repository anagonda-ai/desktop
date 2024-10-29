import os
import gzip
import shutil

def extract_gz(gz_file_path, extract_to_directory):
        
    # Extract the .gz file
    with gzip.open(gz_file_path, 'rb') as f_in:
        output_file_path = os.path.join(extract_to_directory, os.path.splitext(os.path.basename(gz_file_path))[0])
        with open(output_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        print(f"Extracted {gz_file_path} to {output_file_path}")

def process_directory(directory):
    # Iterate through each file in the directory and subdirectories
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".gz"):
                gz_file_path = os.path.join(root, file)
                os.remove(gz_file_path)  # Delete the .gz file
                print(f"Deleted {gz_file_path}")

def main():
    directory = "/groups/itay_mayrose_nosnap/alongonda/full_genomes/phytozome_test/Phytozome/"
    process_directory(directory)

if __name__ == "__main__":
    main()