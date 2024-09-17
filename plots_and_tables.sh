#!/bin/bash

# Directory where the CSV files are located
CSV_DIR="./CSVs"

# List of required files
required_files=(
    "AS_meson_GEVP.csv"
    "AS_meson_GEVP_mix.csv"
    "F_meson_GEVP_mix.csv"
    "F_meson_GEVP.csv"
    "M1_spectral_density_spectrum.csv"
    "M2_spectral_density_spectrum.csv"
    "M3_spectral_density_spectrum.csv"
    "M4_spectral_density_spectrum.csv"
    "M5_spectral_density_spectrum.csv"
)

# Flag to track if all files exist
all_files_exist=true

# Check if each required file exists
for file in "${required_files[@]}"; do
    if [[ ! -f "$CSV_DIR/$file" ]]; then
        echo "Missing required file: $CSV_DIR/$file"
        all_files_exist=false
    fi
done

# Run the commands only if all files exist
if $all_files_exist; then
    python CSVs_to_tables.py
    bash run_plots.sh
else
    echo "One or more required files are missing. Aborting execution."
fi

