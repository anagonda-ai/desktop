import os

def merge_files_unique_ids(base_dir, output_file):
  merged_data = {}
  count_genes = 0
  for dirpath, _, filenames in os.walk(base_dir):
    if "sequences_for_each_port_id.txt" in filenames:
      file_path = os.path.join(dirpath, "sequences_for_each_port_id.txt")

      with open(file_path, "r") as f:
        for line in f:
          if line.startswith('>'):
            protein_id = line.strip()[1:]
            sequence = ""
          else:
            sequence += line.strip()
            if protein_id not in merged_data:
                count_genes += 1
                merged_data[protein_id] = sequence

  with open(os.path.join(base_dir, output_file), "w") as f:
    for protein_id, sequence in merged_data.items():
      f.write(f">{protein_id}\n{sequence}\n")

  print(f"Merged data saved to: {output_file}")
  print(f"{count_genes} unique genes with sequence in plantcyc directory")

if __name__ == "__main__":
  base_dir = "/groups/itay_mayrose_nosnap/alongonda/plantcyc/"
  output_file = "merged_sequences.txt"
  merge_files_unique_ids(base_dir, output_file)