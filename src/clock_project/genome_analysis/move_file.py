import argparse
import os
import shutil
import random
import json

def move_files(file_paths, destination, proportion):
    if not os.path.exists(destination):
        os.makedirs(destination)

    total_files = len(file_paths)
    num_files_to_move = int(total_files * proportion)
    files_to_move = random.sample(file_paths, num_files_to_move)

    for file_path in files_to_move:
        shutil.move(file_path, destination)
        print(f"Moved {file_path} to {destination}")

def main(json_file, destination, proportion):
    with open(json_file, 'r') as f:
        file_paths = json.load(f)
    
    move_files(file_paths, destination, proportion)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Move a proportion of files to a different folder.')
    parser.add_argument('json_file', type=str, help='Path to the JSON file containing the list of file paths to move')
    parser.add_argument('destination', type=str, help='Destination folder to move the files to')
    parser.add_argument('proportion', type=float, help='Proportion of files to move (e.g., 0.2 for 20%)')

    args = parser.parse_args()
    main(args.json_file, args.destination, args.proportion)
