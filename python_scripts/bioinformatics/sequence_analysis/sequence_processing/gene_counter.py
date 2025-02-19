import os

def count_rows_in_files(directory):
    """Counts the number of rows in files named "sequences_for_each_port_id.txt" within a given directory.

    Args:
        directory (str): The path to the directory to search.

    Returns:
        list: A list of tuples containing the subdirectory name and the number of rows in the file.
    """

    with open(output_file, "w") as output:
        for subdir, _, files in os.walk(directory):
            for file in files:
                if file == "sequences_for_each_port_id.txt":
                    file_path = os.path.join(subdir, file)
                    with open(file_path, "r") as f:
                        row_count = len(f.readlines())
                    output.write(f"{subdir.split('/')[-3]}: {row_count / 2}\n")

# Replace with your desired directory path
directory = "/groups/itay_mayrose/alongonda/datasets/plantcyc/all_organisms/"
output_file = "/groups/itay_mayrose/alongonda/datasets/plantcyc/all_organisms_row_count.txt"
# Get the results
count_rows_in_files(directory)