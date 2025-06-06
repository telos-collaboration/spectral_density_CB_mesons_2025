import shutil
import os
import json
import numpy as np

# Optional: Copy relevant directories first
def copy_directory(subdir_name):
    src_dir = os.path.join("..", "input_fit", subdir_name)
    dst_dir = os.path.join("..", "CB_autocorrelation_decay_constant", subdir_name)

    if not os.path.exists(src_dir):
        raise FileNotFoundError(f"Source directory not found: {src_dir}")

    if os.path.exists(dst_dir):
        shutil.rmtree(dst_dir)

    shutil.copytree(src_dir, dst_dir)
    print(f"Copied '{src_dir}' to '{dst_dir}'")

# Step 1: Copy 'metadata' and 'raw_data'
copy_directory("metadata")
copy_directory("raw_data")

# Step 2: Create 'external_data/smeared' directory
base_dir = os.path.join("..", "CB_autocorrelation_decay_constant", "external_data", "smeared")
os.makedirs(base_dir, exist_ok=True)
print(f"Created directory: {base_dir}")

# Step 3: Copy subdirectories from '../JSONs/' to 'smeared'
jsons_src = os.path.join("..", "JSONs")

if not os.path.exists(jsons_src):
    raise FileNotFoundError(f"Source directory not found: {jsons_src}")

for item in os.listdir(jsons_src):
    src_path = os.path.join(jsons_src, item)
    dst_path = os.path.join(base_dir, item)
    if os.path.isdir(src_path):
        shutil.copytree(src_path, dst_path)
        print(f"Copied directory '{src_path}' to '{dst_path}'")

# Step 4: Adjustments for specific (directory, observable) pairs
adjustments = {
    ("Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20", "f_ps"): -0.0004,
    ("Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20", "f_v"): -0.004,
    ("Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20", "f_v"): -0.002,
    ("Sp4b6.5nF2nAS3mF-0.72mAS-1.01T64L32", "f_v"): 0.003,
    ("Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20", "as_ps"): 0.002,
    ("Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20", "as_ps"): 0.002,
    ("Sp4b6.5nF2nAS3mF-0.72mAS-1.01T64L32", "as_ps"): 0.002,
    ("Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20", "as_v"): 0.007,
    ("Sp4b6.5nF2nAS3mF-0.71mAS-1.01T96L20", "as_v"): 0.004,
}

# Bootstrap function
def bootstrap_error(samples, n_resamples=10000):
    resampled_means = [
        np.mean(np.random.choice(samples, size=len(samples), replace=True))
        for _ in range(n_resamples)
    ]
    return np.std(resampled_means)

# Step 5: Apply adjustments and save to modified_ files
for (directory, observable), adjustment in adjustments.items():
    orig_file_path = os.path.join(base_dir, directory, f"meson_extraction_{observable}_samples.json")

    # Load JSON
    with open(orig_file_path, "r") as file:
        data = json.load(file)

    samples_key = f"{observable}_matrix_element"
    if samples_key not in data:
        print(f"Warning: Key '{samples_key}' not found in {orig_file_path}")
        continue

    # Apply adjustment
    adjusted_samples = [x + adjustment for x in data[samples_key]]
    data[samples_key] = adjusted_samples

    # Compute stats
    average = np.mean(adjusted_samples)
    error = bootstrap_error(adjusted_samples)

    # Print summary
    print(f"\nDirectory: {directory}, Observable: {observable}")
    print(f"Adjusted average:        {average:.8f}")
    print(f"Bootstrap error:         {error:.8f}")

    # Save to new file in the same directory
    modified_file_path = os.path.join(base_dir, directory, f"meson_extraction_{observable}_samples.json")
    with open(modified_file_path, "w") as file:
        json.dump(data, file, indent=2)
    print(f"Saved modified file: {modified_file_path}")

