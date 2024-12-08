#!/bin/bash

# Define the directory to search in and the output file
search_dir="/groups/itay_mayrose_nosnap/alongonda/full_genomes/phytozome/Phytozome"  # Update this to your directory
output_file="/groups/itay_mayrose/alongonda/desktop/merged.gff3"  # Update this to your output file

# Clear the output file if it exists
> "$output_file"

# Find all .gff3 files and merge them
find "$search_dir" -type f -name "*.gff3" | while read -r file; do
  # Append the content of each .gff3 file to the output file
  cat "$file" >> "$output_file"
  # Add a newline for separation (optional)
  echo >> "$output_file"
done

echo "All .gff3 files have been merged into $output_file"
