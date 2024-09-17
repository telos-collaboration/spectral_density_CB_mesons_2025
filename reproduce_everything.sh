#!/bin/bash

subdirs=("M1" "M2" "M3" "M4" "M5")
subdirs2=("CSVs")
directory="input_fit"
all_subdirs_present=true



# Check if LaTeX is installed
if ! command -v latex > /dev/null 2>&1; then
    cd CSVs
    echo "LaTeX is not installed. Producing CSVs, no plots nor tables.tex files"
    inner_condition_met=true
    for subdir in "${subdirs2[@]}"; do
        echo "Checking files for pattern: *_meson*"
        count=$(find . -maxdepth 1 -name "*_meson*" | wc -l)
        echo "Found $count subdirectories matching $subdir/*_meson*"
        echo "$count"
        if [ "$count" -ne 5 ]; then
            inner_condition_met=false
            break
        fi
    done
    cd ..

    if [ "$inner_condition_met" = false ]; then
	bash run_plateaus.sh
        bash run_spectral_densities.sh
    else
        bash run_spectral_densities.sh
    fi
    
else
    cd CSVs
    echo "Latex is installed. Running full workflow."
    inner_condition_met=true
    for subdir in "${subdirs2[@]}"; do
        echo "Checking files for pattern: *_meson*"
        count=$(find . -maxdepth 1 -name "*_meson*" | wc -l)
        echo "Found $count subdirectories matching $subdir/*_meson*"
        echo "$count"
        if [ "$count" -ne 5 ]; then
            inner_condition_met=false
            break
        fi
    done
    cd ..

    if [ "$inner_condition_met" = false ]; then
	bash run_plateaus.sh
        bash run_spectral_densities.sh
        python CSVs_to_tables.py
        bash run_plots.sh
    else
        bash run_spectral_densities.sh
        python CSVs_to_tables.py
        bash run_plots.sh
    fi
fi
