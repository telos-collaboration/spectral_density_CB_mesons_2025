#!/bin/bash

# Check if LaTeX is installed
if ! command -v latex > /dev/null 2>&1; then
    source activate
    conda activate snakemake_x86
    bash run_plateaus.sh
    conda activate analysis-env4
    bash run_spectral_densities.sh
    
else
    source activate
    conda activate snakemake_x86
    bash run_plateaus.sh
    conda activate analysis-env4
    bash run_spectral_densities.sh
    CSVs_to_tables_meson_GEVP.py
    CSVs_to_tables_CB_GEVP.py
    CSVs_to_tables_meson_matrixelements.py
    CSVs_to_tables_CB_matrixelements.py
    python renormalise.py
    bash run_plots.sh
fi
