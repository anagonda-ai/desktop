import os
import shutil

def process_fasta_files(base_dir):
    for root, dirs, files in os.walk(base_dir):
        for dir_name in dirs:
            dir_path = os.path.join(root, dir_name)
            if dir_name == "orthology":
                shutil.rmtree(dir_path)
                print(f"Removed directory: {dir_path}")

def main():
    base_dir = "/groups/itay_mayrose/alongonda/full_genomes/phytozome/Phytozome"
    process_fasta_files(base_dir)

if __name__ == "__main__":
    main()