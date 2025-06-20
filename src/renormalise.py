import argparse

import pandas as pd
import re
import os
import numpy as np
import h5py
import math

# Constants
PI2 = math.pi**2

parser = argparse.ArgumentParser()
parser.add_argument("--topology_h5", required=True)
args = parser.parse_args()

# Read the renormalise.csv file
df = pd.read_csv("./metadata/renormalise.csv")

# Compute Z values dynamically using CSV and HDF5 plaquette data
computed_z = {}

with h5py.File(args.topology_h5, "r") as f:
    for _, row in df.iterrows():
        ens = row["Ens"]
        beta = row["beta"]
        C_fund = row["C_fund"]
        C_as = row["C_as"]
        Delta_Sigma1 = row["Delta_Sigma1"]
        Delta_R1 = row["Delta_R1"]
        Delta_R2 = row["Delta_R2"]
        Delta_Lambda = row["Delta_Lambda"]
        Delta_Sigma = row["Delta_Sigma"]

        # Load and average plaquette
        plaq_path = f"{ens}/plaquette"
        if plaq_path not in f:
            raise KeyError(f"Missing {plaq_path} in {args.topology_h5}")
        plaq_values = f[plaq_path][()]
        plaq_avg = np.mean(plaq_values)

        # Compute Zs
        factor_fund = 8 * C_fund / (16 * PI2 * beta * plaq_avg)
        factor_as = 8 * C_as / (16 * PI2 * beta * plaq_avg)

        computed_z[ens] = {
            "Z_PS_fund": 1 + factor_fund * (Delta_Sigma1 + Delta_R1),
            "Z_V_fund": 1 + factor_fund * (Delta_Sigma1 + Delta_R2),
            "Z_A_fund": 1 + factor_fund * (Delta_Sigma1 + Delta_R1),
            "Z_PS_as": 1 + factor_as * (Delta_Sigma1 + Delta_R1),
            "Z_V_as": 1 + factor_as * (Delta_Sigma1 + Delta_R2),
            "Z_A_as": 1 + factor_as * (Delta_Sigma1 + Delta_R1),
            "Z_Lambda": 1
            + (factor_fund / C_fund)
            * ((C_fund + 1 / 2 * C_as) * Delta_Sigma1 + Delta_Lambda),
            "Z_Sigma": 1
            + (factor_fund / C_fund)
            * ((C_fund + 1 / 2 * C_as) * Delta_Sigma1 + Delta_Sigma),
        }

print(computed_z)

# Path to the CSV file
csv_file_path = "./metadata/renormalise.csv"

# Ensembles and file mappings
ensembles = ["M1", "M2", "M3", "M4", "M5"]
tex_files = {
    "M1": "./assets/tables/M1_matrix_meson.tex",
    "M2": "./assets/tables/M2_matrix_meson.tex",
    "M3": "./assets/tables/M3_matrix_meson.tex",
    "M4": "./assets/tables/M4_matrix_meson.tex",
    "M5": "./assets/tables/M5_matrix_meson.tex",
}


# Mapping of table row names to CSV columns
label_to_z_field = {
    "PS": "Z_PS_fund",
    "V": "Z_V_fund",
    "AV": "Z_A_fund",
    "ps": "Z_PS_as",
    "v": "Z_V_as",
    "av": "Z_A_as",
}


# Loop over each ensemble and associated file
for ens, tex_file in tex_files.items():
    # Get the row for the current ensemble
    z_values = computed_z[ens]
    # Read the .tex file line by line
    with open(tex_file, "r") as file:
        lines = file.readlines()

    # Process and update each line with the corresponding Z value
    updated_lines = []

    for line in lines:
        line_modified = False  # Track if the line was modified

        # Check if the line contains a label we want to replace
        for label, z_field in label_to_z_field.items():
            # Retrieve the Z value for this label
            z_value = z_values[z_field]
            print(f"Processing label: {label}, Z value: {z_value}")

            # Check if the line contains the label (e.g., "PS", "V", etc.)
            if re.search(rf"^\s*{label}\s*&", line):
                # Find all matches of numbers with uncertainty in the line
                matches = re.finditer(r"([-+]?\d*\.\d+)\((\d+)\)", line)

                updated_line = line
                for match in matches:
                    # Extract main value and uncertainty
                    main_value = float(match.group(1))
                    uncertainty = match.group(2)

                    # Determine the number of decimal places in the original value
                    original_value_str = match.group(
                        1
                    )  # This is the original string of the main value
                    decimal_places = (
                        len(original_value_str.split(".")[1])
                        if "." in original_value_str
                        else 0
                    )

                    # Perform the multiplication on the main value
                    new_main_value = main_value * z_value

                    # Dynamically format the result based on the original number of decimal places
                    format_str = f"{{:.{decimal_places}f}}({uncertainty})"
                    new_value_str = format_str.format(new_main_value)

                    # Replace the original value with the new calculated value
                    updated_line = updated_line.replace(
                        f"{main_value:.{decimal_places}f}({uncertainty})", new_value_str
                    )

                    print(
                        f"Replaced '{main_value:.{decimal_places}f}({uncertainty})' with '{new_value_str}' in line for label '{label}'"
                    )

                # Append the updated line to the list and mark it as modified
                updated_lines.append(updated_line)
                line_modified = True
                break  # Stop after finding and modifying the line for this label

        if not line_modified:
            # If no label matches, keep the line as is
            updated_lines.append(line)

    # Save the updated lines back to a new .tex file
    output_file = f"./assets/tables/renormalised_{ens}_matrix_meson.tex"
    with open(output_file, "w") as file:
        file.writelines(updated_lines)

    print(f"Updated {tex_file} and saved as {output_file}")


ensembles = ["M1", "M2", "M3", "M4", "M5"]

for ens in ensembles:
    z_values = computed_z[ens]
    cb_file = f"./assets/tables/{ens}_matrix_CB.tex"
    cb_output = f"./assets/tables/renormalised_{ens}_matrix_CB.tex"

    if not os.path.exists(cb_file):
        print(f"Missing: {cb_file}")
        continue

    with open(cb_file, "r") as file:
        lines = file.readlines()

    updated_lines = []

    for line in lines:
        original_line = line
        # Only modify lines containing data (skip headers, hline etc.)
        if re.match(r"^\s*\$", line.strip()):
            # Decide Z factor
            if "Lambda" in line:
                z = z_values["Z_Lambda"]
            elif "Sigma" in line:
                z = z_values["Z_Sigma"]
            else:
                updated_lines.append(line)
                continue

            # Match numbers with uncertainties: only in columns 2-4
            matches = list(re.finditer(r"([-+]?\d*\.\d+)\((\d+)\)", line))
            if len(matches) < 3:
                updated_lines.append(line)
                continue

            updated_line = line
            for i in range(3):  # Only apply to the first three matches
                match = matches[i]
                main_val = float(match.group(1))
                unc = match.group(2)
                dec = len(match.group(1).split(".")[1])
                new_val = main_val * z
                formatted = f"{new_val:.{dec}f}({unc})"
                original = match.group(0)
                updated_line = updated_line.replace(original, formatted, 1)

            updated_lines.append(updated_line)
        else:
            updated_lines.append(line)

    with open(cb_output, "w") as f:
        f.writelines(updated_lines)

    print(f"Renormalised {cb_file} and saved to {cb_output}")
