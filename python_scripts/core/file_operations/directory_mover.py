import os
import shutil

def move_directories(base_dir, new_dir, exclude_dir):
    # Create the new directory if it doesn't exist
    os.makedirs(new_dir, exist_ok=True)
    
    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        if os.path.isdir(item_path) and item != exclude_dir:
            shutil.move(item_path, new_dir)
            print(f"Moved directory: {item_path} to {new_dir}")

def main():
    base_dir = "/groups/itay_mayrose/alongonda/datasets/full_genomes/plaza"
    new_dir = "/groups/itay_mayrose/alongonda/datasets/full_genomes/plaza/organisms"
    exclude_dir = "processed_annotations"
    
    move_directories(base_dir, new_dir, exclude_dir)

if __name__ == "__main__":
    main()