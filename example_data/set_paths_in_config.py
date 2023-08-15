#!/usr/bin/env python3
import argparse
import os


def modify_config_file(data_dir, config_file):
    # Ensure the path is absolute
    if not os.path.isabs(data_dir):
        raise ValueError("The provided data path is not an absolute path.")

    # Ensure data_dir ends with a '/'
    if not data_dir.endswith('/'):
        data_dir += '/'

    # Read the configuration file
    with open(config_file, "r") as file:
        lines = file.readlines()

    # Modify the relevant lines
    for i, line in enumerate(lines):
        if line.startswith("data_dir ="):
            lines[i] = f"data_dir = {data_dir}\n"
        elif line.startswith("segment_split ="):
            filename = os.path.basename(line.split('=')[1].strip())
            lines[i] = f"segment_split = {os.path.join(data_dir, filename)}\n"
        elif line.startswith("centromeres ="):
            filename = os.path.basename(line.split('=')[1].strip())
            lines[i] = f"centromeres = {os.path.join(data_dir, filename)}\n"

    # Save the modified file
    with open(config_file, "w") as file:
        file.writelines(lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Modifies the provided stg_ini file to set paths based on the provided data folder.')
    parser.add_argument('-d', '--data_dir', required=True, help='Absolute path to the data folder.')
    parser.add_argument('-c', '--file', required=True, help='Path to the stg_ini configuration file.')
    args = parser.parse_args()

    modify_config_file(args.data_dir, args.file)
