#!/bin/bash

set -eu

# Check if LaTeX is installed
if ! command -v latex > /dev/null 2>&1; then
    bash run_plateaus.sh
    cd ./CB_autocorrelation_decay_constant
    bash main.sh
    cd ..
    bash run_spectral_densities.sh
    python CSVs_to_tables_meson_GEVP.py
    python CSVs_to_tables_CB_GEVP.py
    python CSVs_to_tables_meson_matrixelements.py
    python CSVs_to_tables_CB_matrixelements.py
    python renormalise.py
else
    bash run_plateaus.sh
    cd ./CB_autocorrelation_decay_constant
    bash main.sh
    cd ..
    bash run_spectral_densities.sh
    python CSVs_to_tables_meson_GEVP.py
    python CSVs_to_tables_CB_GEVP.py
    python CSVs_to_tables_meson_matrixelements.py
    python CSVs_to_tables_CB_matrixelements.py
    python renormalise.py
    bash run_plots.sh
fi
