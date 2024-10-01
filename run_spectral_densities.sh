#!/bin/bash

bash check_latex.sh

directory="input_fit"
subdirs=("M1" "M2" "M3" "M4" "M5")
all_subdirs_present=true
min_size_mb=6
min_size_kb=150

for subdir in "${subdirs[@]}"; do
    if [ ! -d "$directory/$subdir" ]; then
        all_subdirs_present=false
        break
    fi
    # Check if the directory size is greater than 6MB
    dir_size=$(du -sm "$directory/$subdir" | cut -f1)
    if [[ -z "$dir_size" || ! "$dir_size" =~ ^[0-9]+$ ]]; then
        echo "Error: Unable to get size for $directory/$subdir"
        all_subdirs_present=false
        break
    fi
    if [ "$dir_size" -le "$min_size_mb" ]; then
        all_subdirs_present=false
        break
    fi
done

cd lsd_out

if [ "$all_subdirs_present" = false ]; then
    inner_condition_met=true
    for subdir in "${subdirs[@]}"; do
        echo "Checking subdirectories for pattern: ${subdir}_*"
        count=$(find . -maxdepth 1 -type d -name "${subdir}_*" | wc -l)
        echo "Found $count subdirectories matching ${subdir}_*"
        
        if [ "$count" -ne 72 ]; then
            inner_condition_met=false
            break
        fi

        # Check if each matching subdirectory has at least 150KB
        for match in $(find . -maxdepth 1 -type d -name "${subdir}_*"); do
            subdir_size=$(du -sk "$match" | cut -f1)
            if [ "$subdir_size" -lt "$min_size_kb" ]; then
                echo "Subdirectory $match is smaller than ${min_size_kb}KB"
                inner_condition_met=false
                break
            fi
        done

        if [ "$inner_condition_met" = false ]; then
            break
        fi
    done

    if [ "$inner_condition_met" = false ]; then
        python analyse_data_mesons.py
        python analyse_data_CB.py
        python print_samples_mesons.py
	python print_samples_CB.py
        python fit_data_mesons.py
        python fit_data_CB.py
        python simultaneous_fits_mesons.py
        python simultaneous_fits_CB.py
    else
        python print_samples_mesons.py
        python print_samples_CB.py
        python fit_data_mesons.py
        python fit_data_CB.py
        python simultaneous_fits_mesons.py
        python simultaneous_fits_CB.py
    fi

else
    python fit_data_mesons.py
    python fit_data_CB.py
    python simultaneous_fits_mesons.py
    python simultaneous_fits_CB.py
fi

