import shutil
import os
import json
import numpy as np
import csv
import re

def copy_directory(subdir_name):
    src_dir = os.path.join("..", "input_fit", subdir_name)
    dst_dir = os.path.join("..", "CB_autocorrelation_decay_constant", subdir_name)

    if not os.path.exists(src_dir):
        raise FileNotFoundError(f"Source directory not found: {src_dir}")

    if os.path.exists(dst_dir):
        shutil.rmtree(dst_dir)

    shutil.copytree(src_dir, dst_dir)
    print(f"Copied '{src_dir}' to '{dst_dir}'")

copy_directory("metadata")
copy_directory("raw_data")

base_dir = os.path.join("..", "CB_autocorrelation_decay_constant", "external_data", "smeared")
os.makedirs(base_dir, exist_ok=True)
print(f"Created directory: {base_dir}")

jsons_src = os.path.join("..", "JSONs")

if not os.path.exists(jsons_src):
    raise FileNotFoundError(f"Source directory not found: {jsons_src}")

for item in os.listdir(jsons_src):
    src_path = os.path.join(jsons_src, item)
    dst_path = os.path.join(base_dir, item)
    if os.path.isdir(src_path):
        shutil.copytree(src_path, dst_path)
        print(f"Copied directory '{src_path}' to '{dst_path}'")




# Load metadata
plateau_data = {}
with open('../input_fit/metadata_plateaus.csv', newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        key = (row['ensemble_name'].strip(), row['channel'].strip())
        plateau_data[key] = (int(row['plateau_start']), int(row['plateau_end']))

# Read and update main.sh
with open('../CB_autocorrelation_decay_constant/main.sh', 'r') as file:
    lines = file.readlines()

updated_lines = []
pattern = re.compile(
    r'^(python3\s+src_py/mass_wall\.py\s+--ensemble_name\s+)(\S+)\s+--plateau_start\s+\d+\s+--plateau_end\s+\d+(.*--channel\s+)(\S+)(\s+data_assets/wall_correlators\.hdf5.*)$'
)

for line in lines:
    match = pattern.match(line)
    if match:
        pre, ensemble_name, mid, channel, post = match.groups()
        key = (ensemble_name, channel)
        if key in plateau_data:
            start, end = plateau_data[key]
            new_line = f"{pre}{ensemble_name} --plateau_start {start} --plateau_end {end}{mid}{channel}{post}\n"
            updated_lines.append(new_line)
        else:
            print(f"Warning: No metadata found for {key}, keeping original line.")
            updated_lines.append(line)
    else:
        updated_lines.append(line)


with open('../CB_autocorrelation_decay_constant/main.sh', 'w') as file:
    file.writelines(updated_lines)

