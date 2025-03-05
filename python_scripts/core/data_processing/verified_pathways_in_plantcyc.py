import os
import csv
import re

# Define base directory containing all organism databases
BASE_DIR = "/groups/itay_mayrose/alongonda/datasets/plantcyc/all_organisms"
OUTPUT_DIR = os.path.join(BASE_DIR, "output")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Regex pattern to detect MetaCyc DBLINKS
METACYC_PATTERN = re.compile(r'DBLINKS\s+-\s+\(METACYC\s+"([^"]+)"\)')

def extract_metacyc_pathways(file_path):
    """ Parses pathways.dat and extracts pathways with MetaCyc DBLINKS. """
    metacyc_pathways = []

    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        pathway_id = None
        common_name = None
        has_metacyc_dblink = False
        metacyc_id = None

        for line in f:
            line = line.strip()

            # Start of a new pathway entry
            if line.startswith("UNIQUE-ID"):
                # Save the previous pathway if it had a MetaCyc DBLINK
                if pathway_id and has_metacyc_dblink:
                    metacyc_pathways.append([pathway_id, common_name or "Unknown", metacyc_id])

                # Reset for new pathway
                pathway_id = line.split(" - ")[1] if " - " in line else None
                common_name = None
                has_metacyc_dblink = False
                metacyc_id = None

            # Extract pathway name
            elif line.startswith("COMMON-NAME"):
                common_name = line.split(" - ", 1)[1] if " - " in line else None

            # Check for MetaCyc DBLINKS
            elif line.startswith("DBLINKS"):
                match = METACYC_PATTERN.search(line)
                if match:
                    has_metacyc_dblink = True
                    metacyc_id = match.group(1)

        # Add last pathway if it has a MetaCyc DBLINK
        if pathway_id and has_metacyc_dblink:
            metacyc_pathways.append([pathway_id, common_name or "Unknown", metacyc_id])

    return metacyc_pathways

num_of_pathways = 0
# Walk through all directories to find pathways.dat
for root, dirs, files in os.walk(BASE_DIR):
    for file in files:
        if file == "pathways.dat":
            pathways_file = os.path.join(root, file)

            # Extract correct organism name (two levels up from `data/`)
            organism_name = os.path.basename(os.path.dirname(os.path.dirname(root)))

            print(f"Processing {organism_name}...")

            # Extract ONLY pathways linked to MetaCyc
            metacyc_pathways = extract_metacyc_pathways(pathways_file)

            if metacyc_pathways:  # Only save if there are MetaCyc pathways
                # Define output file for this organism
                output_file = os.path.join(OUTPUT_DIR, f"{organism_name}_metacyc_pathways.csv")

                # Write results to CSV
                with open(output_file, "w", newline="", encoding="utf-8") as f:
                    writer = csv.writer(f)
                    writer.writerow(["Pathway ID", "Common Name", "MetaCyc ID"])
                    writer.writerows(metacyc_pathways)

                print(f"✅ Saved {len(metacyc_pathways)} MetaCyc pathways to {output_file}")
                num_of_pathways += len(metacyc_pathways)
            else:
                print(f"❌ No MetaCyc pathways found for {organism_name}")

print(f"📊 Total MetaCyc pathways found: {num_of_pathways}")
print("🎉 Done processing all databases!")
