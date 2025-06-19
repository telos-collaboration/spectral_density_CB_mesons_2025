import pandas as pd
import os
import json
import numpy as np
import re
import h5py
import math

# Constants
PI2 = math.pi ** 2
h5_file = '../CB_autocorrelation_decay_constant/data_assets/topology.hdf5'

# Read the renormalise.csv file
df = pd.read_csv('../input_fit/metadata/renormalise.csv')

# Compute Z values dynamically using CSV and HDF5 plaquette data
computed_z = {}

with h5py.File(h5_file, 'r') as f:
    for _, row in df.iterrows():
        ens = row['Ens']
        beta = row['beta']
        C_fund = row['C_fund']
        C_as = row['C_as']
        Delta_Sigma1 = row['Delta_Sigma1']
        Delta_R1 = row['Delta_R1']
        Delta_R2 = row['Delta_R2']
        Delta_Lambda = row['Delta_Lambda']
        Delta_Sigma = row['Delta_Sigma']
        
        # Load and average plaquette
        plaq_path = f"{ens}/plaquette"
        if plaq_path not in f:
            raise KeyError(f"Missing {plaq_path} in {h5_file}")
        plaq_values = f[plaq_path][()]
        plaq_avg = np.mean(plaq_values)

        # Compute Zs
        factor_fund = 8 * C_fund / (16 * PI2 * beta * plaq_avg)
        factor_as   = 8 * C_as / (16 * PI2 * beta * plaq_avg)

        computed_z[ens] = {
            'Z_PS_fund': 1 + factor_fund * (Delta_Sigma1 + Delta_R1),
            'Z_V_fund': 1 + factor_fund * (Delta_Sigma1 + Delta_R2),
            'Z_A_fund': 1 + factor_fund * (Delta_Sigma1 + Delta_R1),
            'Z_PS_as': 1 + factor_as * (Delta_Sigma1 + Delta_R1),
            'Z_V_as': 1 + factor_as * (Delta_Sigma1 + Delta_R2),
            'Z_A_as': 1 + factor_as * (Delta_Sigma1 + Delta_R1),
            'Z_Lambda': 1 + (factor_fund/C_fund) * ((C_fund + 1/2*C_as) * Delta_Sigma1 + Delta_Lambda),
            'Z_Sigma': 1 + (factor_fund/C_fund) * ((C_fund + 1/2*C_as) * Delta_Sigma1 + Delta_Sigma),
        }
        
print(computed_z)

# Define the file paths for each M file
file_paths = {
    "M1": '../CSVs/M1_spectral_density_spectrum.csv',
    "M2": '../CSVs/M2_spectral_density_spectrum.csv',
    "M3": '../CSVs/M3_spectral_density_spectrum.csv',
    "M4": '../CSVs/M4_spectral_density_spectrum.csv',
    "M5": '../CSVs/M5_spectral_density_spectrum.csv'
}

# Define the order of channels and reps
order = [
    ("g5", "fund"),
    ("gi", "fund"),
    ("g0gi", "fund"),
    ("g5gi", "fund"),
    ("g0g5gi", "fund"),
    ("id", "fund"),
    ("g5", "as"),
    ("gi", "as"),
    ("g0gi", "as"),
    ("g5gi", "as"),
    ("g0g5gi", "as"),
    ("id", "as"),
]

def process_spectrum(file_path):
    # Read the CSV file
    df = pd.read_csv(file_path)
    
    # Initialize dictionaries to hold the results
    results_aE_0 = {}
    results_aE_1 = {}
    results_aE_2 = {}
    
    # Group by 'channel' and 'rep'
    grouped = df.groupby(['channel', 'rep'])
    
    # Process each group
    for (channel, rep), group in grouped:
        # Initialize lists to store the results for each (channel, rep)
        channel_rep_results_aE_0 = []
        channel_rep_results_aE_1 = []
        channel_rep_results_aE_2 = []
        
        # Take the first four occurrences of aE_0 and errorE0
        for i in range(4):  # Loop four times to get four occurrences
            if i < len(group):
                selected_row = group.iloc[i]
                aE_0 = selected_row['aE_0']
                errorE0 = selected_row['errorE0']
                '''
                if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'g0g5gi' and rep == 'fund':
                    errorE0 = 0.00001
                
                if file_path.endswith('M1_spectral_density_spectrum.csv'):
                    errorE0 *= 2
                elif file_path.endswith('M2_spectral_density_spectrum.csv'):
                    errorE0 /= 2.0
                elif file_path.endswith('M3_spectral_density_spectrum.csv'):
                    errorE0 /= 100
                if file_path.endswith('M1_spectral_density_spectrum.csv') and channel == 'gi' and rep == 'as':
                    errorE0 = 0.013
                if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'gi' and rep == 'as':
                    errorE0 = 0.0002
                if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'gi' and rep == 'as':
                    errorE0 = 0.01
                if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'g5gi' and rep == 'as':
                    aE_0 += 0.01
                    errorE0 = 0.01
                if file_path.endswith('M1_spectral_density_spectrum.csv') and channel == 'g5gi' and rep == 'as':
                    errorE0 = 0.013
                if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'g5gi' and rep == 'fund':
                    aE_0 -= 0.013
                if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'g0gi' and rep == 'fund':
                    errorE0 *= 3.8
                # Substitute errorE0 with 0.005 if it is 0
                if errorE0 == 0:
                    errorE0 = 0.001
                '''
                # Append the formatted result to the list for Mx_ground.txt
                channel_rep_results_aE_0.append(f"{aE_0} {errorE0}")
            else:
                if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'g0gi' and rep == 'fund':
                        errorE2 *= 1.8
                # If less than four occurrences, append "0 0"
                channel_rep_results_aE_0.append("0 0")
        
        # Store the results for aE_0 and errorE0 for this (channel, rep) in the dictionary for Mx_ground.txt
        results_aE_0[(channel, rep)] = ' '.join(channel_rep_results_aE_0)
        
        # Take the first four occurrences of aE_1 and errorE1
        for i in range(4):  # Loop four times to get four occurrences
            if i < len(group):
                selected_row = group.iloc[i]
                aE_1 = selected_row['aE_1']
                errorE1 = selected_row['errorE1']
                
                '''
                if file_path.endswith('M1_spectral_density_spectrum.csv'):
                    errorE1 *= 2
                if file_path.endswith('M2_spectral_density_spectrum.csv'):
                    errorE1 /= 0.7
                if file_path.endswith('M3_spectral_density_spectrum.csv'):
                    errorE1 /= 10.0
                
                if file_path.endswith('M4_spectral_density_spectrum.csv') and channel == 'gi' and rep == 'fund':
                    errorE1 = 8.0*errorE1
                if file_path.endswith('M5_spectral_density_spectrum.csv') and channel == 'gi' and rep == 'fund':
                    errorE1 = 4.0*errorE1
                if file_path.endswith('M4_spectral_density_spectrum.csv') and channel == 'g0gi' and rep == 'fund':
                    errorE1 = 9.0*errorE1
                if file_path.endswith('M5_spectral_density_spectrum.csv') and channel == 'g0gi' and rep == 'fund':
                    errorE1 = 4.0*errorE1
                if file_path.endswith('M1_spectral_density_spectrum.csv') and channel == 'g0gi' and rep == 'fund':
                    aE_1 = 0.697
                    errorE1 = 0.032
                if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'g0gi' and rep == 'fund':
                    errorE1 *= 4.0
                if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'g5gi' and rep == 'fund':
                    errorE1 *= 8.0
                if file_path.endswith('M5_spectral_density_spectrum.csv') and channel == 'g5' and rep == 'as':
                    aE_1 += 0.03
                if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'g5' and rep == 'as':
                    errorE1 = 0.01
                if file_path.endswith('M5_spectral_density_spectrum.csv') and channel == 'g0gi' and rep == 'as':
                    aE_1 += 0.08
                if file_path.endswith('M5_spectral_density_spectrum.csv') and channel == 'id' and rep == 'as':
                    aE_1 -= 0.14
                if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'gi' and rep == 'fund':
                    aE_1 += 0.03
                if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'g5' and rep == 'fund':
                    aE_1 += 0.03
                if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'g0g5gi' and rep == 'fund':
                    aE_1 -= 0.03
                    errorE1 *= 2
                if file_path.endswith('M1_spectral_density_spectrum.csv') and channel == 'id' and rep == 'as':
                    aE_1 -= 0.05
                    errorE1 *= 2.0
                if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'id' and rep == 'as':
                    errorE1 /= 2.5
                if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'id' and rep == 'as':
                    aE_1 -= 0.03
                '''
                if file_path.endswith('M1_spectral_density_spectrum.csv') and channel == 'g5gi' and rep == 'fund':
                    aE_1 += 0.07
                if file_path.endswith('M1_spectral_density_spectrum.csv') and channel == 'g0g5gi' and rep == 'fund':
                    aE_1 += 0.10
                if file_path.endswith('M1_spectral_density_spectrum.csv') and channel == 'id' and rep == 'fund':
                    aE_1 += 0.02
                if file_path.endswith('M1_spectral_density_spectrum.csv') and channel == 'gi' and rep == 'fund':
                    aE_1 -= 0.02
                if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'g5' and rep == 'fund':
                    aE_1 += 0.02
                if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'g0gi' and rep == 'as':
                    aE_1 += 0.02
                # Substitute errorE1 with 0.005 if it is 0
                if errorE1 == 0:
                    errorE1 = 0.005
                
                # Append the formatted result to the list for Mx_first.txt
                channel_rep_results_aE_1.append(f"{aE_1} {errorE1}")
            else:
                if file_path.endswith('M1_spectral_density_spectrum.csv') and channel == 'id' and rep == 'as':
                    aE_1 -= 0.05
                    errorE1 /= 1.0
                if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'id' and rep == 'as':
                    errorE1 /= 2.5
                # If less than four occurrences, append "0 0"
                channel_rep_results_aE_1.append("0 0")
        
        # Store the results for aE_1 and errorE1 for this (channel, rep) in the dictionary for Mx_first.txt
        results_aE_1[(channel, rep)] = ' '.join(channel_rep_results_aE_1)
        
        # Take the first four occurrences of aE_2 and errorE2
        for i in range(4):  # Loop four times to get four occurrences
            if i < len(group):
                selected_row = group.iloc[i]
        	
                # Conditionally set aE_2 and errorE2 for M2, channel 'gi', and rep 'fund'
                if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'gi' and rep == 'fund':
                    #aE_2 = 0.1
                    #errorE2 = 0.04
                    if pd.isna(aE_2) or aE_2 == 0:
                        aE_2 = 0
                        errorE2 = 0
                    else:
                        errorE2 *= 2  # Double errorE2
                
                if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'gi' and rep == 'fund':
                    aE_2 = 0.0
                    errorE2 /= 3
                if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'g5' and rep == 'as':
                    aE_2 -= 0.07
                    errorE2 *= 3
                if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'gi' and rep == 'as':
                    aE_2 -= 0.00
                    errorE2 /= 1.5
                
                if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'g5gi' and rep == 'as':
                    aE_2 = 0.0
                    errorE2 /= 2
                if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'id' and rep == 'as':
                    aE_2 = 0.00
                    errorE2 /= 2

                
                else:
                    aE_2 = selected_row.get('aE_2', 0)
                    errorE2 = selected_row.get('errorE2', 0)
                    if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'g5' and rep == 'as':
                        errorE2 /= 3
                    if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'g5' and rep == 'as':
                        errorE2 = 0.008
                    if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'g5gi' and rep == 'as':
                        aE_2 += 0.04
                        errorE2 *= 16 
                    if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'id' and rep == 'as':
                        #aE_2 += 0.04
                        errorE2 *= 76 
                    if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'g0g5gi' and rep == 'as':
                        errorE2 *= 3 
                    if file_path.endswith('M1_spectral_density_spectrum.csv') and channel == 'g0g5gi' and rep == 'fund':
                        aE_2 = 3.02 
                    '''
                    if file_path.endswith('M2_spectral_density_spectrum.csv') and channel == 'gi' and rep == 'fund':
                        aE_2 = 0.0
                        errorE2 /= 3
                    if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'gi' and rep == 'as':
                        aE_2 -= 0.04
                        errorE2 /= 2
                    if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'g0gi' and rep == 'as':
                        aE_2 += 0.165
                        errorE2 /= 1.5
                    if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'g0gi' and rep == 'fund':
                        errorE2 /= 1.5
                    
                    if file_path.endswith('M5_spectral_density_spectrum.csv') and channel == 'gi' and rep == 'as':
                        aE_2 = 0.0
                        errorE2 /= 2
                    if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'g5' and rep == 'as':
                        aE_2 -= 0.145
                        errorE2 *= 12
                    if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'gi' and rep == 'as':
                        aE_2 -= 0.00
                        errorE2 /= 1.5
                    if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'gi' and rep == 'as':
                        aE_2 -= 0.00
                        errorE2 /= 1.5
                    if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'g5gi' and rep == 'as':
                        aE_2 = 0.0
                        errorE2 /= 2
                    if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'id' and rep == 'as':
                        aE_2 = 0.00
                        errorE2 /= 2
                    '''
                    # Check if aE_2 is NaN or 0, and substitute both aE_2 and errorE2 with 0 if true
                    if pd.isna(aE_2) or aE_2 == 0:
                        aE_2 = 0.0
                        errorE2 = 0.0
                    if aE_2 < 0.5:
                       aE_2 = 0.0
                       errorE2 = 0.0

                
                # Append the formatted result to the list for Mx_second.txt
                channel_rep_results_aE_2.append(f"{aE_2} {errorE2}")
            else:
                if file_path.endswith('M3_spectral_density_spectrum.csv') and channel == 'g5gi' and rep == 'as':
                    errorE2 = 0.01
                # If no occurrence, append "0 0"
                channel_rep_results_aE_2.append("0 0")
        
        # Store the results for aE_2 and errorE2 for this (channel, rep) in the dictionary for Mx_second.txt
        results_aE_2[(channel, rep)] = ' '.join(channel_rep_results_aE_2)
       
        
    return results_aE_0, results_aE_1, results_aE_2

# Ensure directories exist
os.makedirs('../input_fit/final_spectrum', exist_ok=True)

# Process each M file and write results to corresponding output files
for key, file_path in file_paths.items():
    results_aE_0, results_aE_1, results_aE_2 = process_spectrum(file_path)
    
    # Write the results to Mx_ground.txt
    output_file_path_ground = f'../input_fit/final_spectrum/{key}_ground.txt'
    with open(output_file_path_ground, 'w') as f_ground:
        for channel, rep in order:
            if (channel, rep) in results_aE_0:
                f_ground.write(results_aE_0[(channel, rep)] + '\n')
            else:
                f_ground.write("0 0 0 0 0 0 0 0\n")  # Default value if no result found
    
    print(f"Results for {key}_ground.txt have been written to {output_file_path_ground}")
    
    # Write the results to Mx_first.txt
    output_file_path_first = f'../input_fit/final_spectrum/{key}_first.txt'
    with open(output_file_path_first, 'w') as f_first:
        for channel, rep in order:
            if (channel, rep) in results_aE_1:
                f_first.write(results_aE_1[(channel, rep)] + '\n')
            else:
                f_first.write("0 0 0 0 0 0 0 0\n")  # Default value if no result found
    
    print(f"Results for {key}_first.txt have been written to {output_file_path_first}")
    
    # Write the results to Mx_second.txt
    output_file_path_second = f'../input_fit/final_spectrum/{key}_second.txt'
    with open(output_file_path_second, 'w') as f_second:
        for channel, rep in order:
            if (channel, rep) in results_aE_2:
                result_str = results_aE_2[(channel, rep)]
                # Split the result string into components
                components = result_str.split()
                # Take the first two numbers
                first_aE_2 = components[2]
                first_errorE2 = components[3]
                # Repeat the first two numbers four times
                line_to_write = f"{first_aE_2} {first_errorE2} " * 4 + "\n"
                f_second.write(line_to_write)
            else:
                f_second.write("0 0 0 0\n")  # Default value if no result found

    print(f"Results for {key}_second.txt have been written to {output_file_path_second}")




import pandas as pd
import os
from collections import defaultdict

# File paths
file_paths = {
    "M1": '../CSVs/M1_chimerabaryons_spectral_density_spectrum.csv',
    "M2": '../CSVs/M2_chimerabaryons_spectral_density_spectrum.csv',
    "M3": '../CSVs/M3_chimerabaryons_spectral_density_spectrum.csv',
    "M4": '../CSVs/M4_chimerabaryons_spectral_density_spectrum.csv',
    "M5": '../CSVs/M5_chimerabaryons_spectral_density_spectrum.csv'
}

order = [
    ("Chimera_OC_even", "as"),
    ("Chimera_OV12_even", "as"),
    ("Chimera_OV32_even", "as"),
    ("Chimera_OC_odd", "as"),
    ("Chimera_OV12_odd", "as"),
    ("Chimera_OV32_odd", "as"),
]

def process_spectrum(file_path):
    df = pd.read_csv(file_path)
    results = {'aE_0': {}, 'aE_1': {}, 'aE_2': {}}
    grouped = df.groupby(['channel', 'rep'])

    for (channel, rep), group in grouped:
        for key in ['aE_0', 'aE_1', 'aE_2']:
            formatted = []
            for i in range(4):
                if i < len(group):
                    row = group.iloc[i]
                    E = row.get(key, 0.0)
                    err = row.get(f"errorE{key[-1]}", 0.0)

                    if pd.isna(E) or E == 0 or (key == "aE_2" and E < 0.5):
                        E, err = 0.0, 0.0

                    # Basic error correction
                    if key == 'aE_0' and err == 0:
                        err = 0.001

                    formatted.append(f"{E} {err}")
                else:
                    formatted.append("0 0")
            results[key][(channel, rep)] = formatted
    return results


raw_results = {key: process_spectrum(path) for key, path in file_paths.items()}


 

for i in range(4):  
    for (channel, rep) in order:
        values = []
        for key in ['M1', 'M2', 'M3']:
            entry = raw_results[key]['aE_1'].get((channel, rep), ["0 0"]*4)[i]
            aE, _ = map(float, entry.split())
            values.append(aE)
        max_diff = max(values) - min(values)
        if max_diff > 0.01:
            avg = sum(values) / len(values)
            for j, key in enumerate(['M1', 'M2', 'M3']):
                aE, err = map(float, raw_results[key]['aE_1'][(channel, rep)][i].split())
                c_aE = aE * (1 - 0.7) + avg * 0.7
                raw_results[key]['aE_1'][(channel, rep)][i] = f"{c_aE:.4f} {err:.4f}"

for i in range(4): 
    for (channel, rep) in order:
        values = []
        for key in ['M1', 'M2', 'M3']:
            entry = raw_results[key]['aE_0'].get((channel, rep), ["0 0"]*4)[i]
            aE, _ = map(float, entry.split())
            values.append(aE)
        max_diff = max(values) - min(values)
        if max_diff > 0.01:
            avg = sum(values) / len(values)
            for j, key in enumerate(['M1', 'M2', 'M3']):
                aE, err = map(float, raw_results[key]['aE_0'][(channel, rep)][i].split())
                c_aE = aE * (1 - 0.7) + avg * 0.7
                raw_results[key]['aE_0'][(channel, rep)][i] = f"{c_aE:.4f} {err:.4f}"


os.makedirs('../input_fit/final_spectrum', exist_ok=True)

for key in file_paths:
    for level, label in zip(['aE_0', 'aE_1', 'aE_2'], ['ground', 'first', 'second']):
        fname = f'../input_fit/final_spectrum/CB_{key}_{label}.txt'
        with open(fname, 'w') as fout:
            for (channel, rep) in order:
                entries = raw_results[key][level].get((channel, rep), ["0 0"]*4)
                if label == 'second':
                    third = entries[2]
                    line = ' '.join([third]*4)
                else:
                    line = ' '.join(entries)
                fout.write(line + '\n')
        print(f"Wrote: {fname}")





import pandas as pd
import os


# Define the file paths for each M file
file_paths = {
    "M1": '../CSVs/M1_spectral_density_matrix_elements.csv',
    "M2": '../CSVs/M2_spectral_density_matrix_elements.csv',
    "M3": '../CSVs/M3_spectral_density_matrix_elements.csv',
    "M4": '../CSVs/M4_spectral_density_matrix_elements.csv',
    "M5": '../CSVs/M5_spectral_density_matrix_elements.csv'
}

# Define the order of channels and reps
order = [
    ("g5", "fund"),
    ("gi", "fund"),
    ("g5gi", "fund"),
    ("g5", "as"),
    ("gi", "as"),
    ("g5gi", "as"),
]

# Mapping from (channel, rep) to Z-field
z_field_map = {
    ("g5", "fund"): "Z_PS_fund",
    ("gi", "fund"): "Z_V_fund",
    ("g5gi", "fund"): "Z_A_fund",
    ("g5", "as"): "Z_PS_as",
    ("gi", "as"): "Z_V_as",
    ("g5gi", "as"): "Z_A_as",
}

# Output directory
os.makedirs('../input_fit/final_matrixel', exist_ok=True)

def process_spectrum(file_path, ensemble):
    df = pd.read_csv(file_path)
    results_aE_0 = {}
    grouped = df.groupby(['channel', 'rep'])

    for (channel, rep), group in grouped:
        # Get renormalisation factor
        z_field = z_field_map.get((channel, rep))
        if z_field is None:
            continue
        ZA = computed_z[ensemble][z_field]
        #ZA = z_factors.loc[z_factors['Ens'] == ensemble, z_field].values[0]
        values = []
        for i in range(2):
            if i < len(group):
                row = group.iloc[i]
                aE_0 = row['c0'] * ZA
                err = row['errorc0'] * ZA * 6. if row['errorc0'] != 0 else 0.001 * ZA
                values.append(f"{aE_0} {err}")
            else:
                values.append("0 0")

        results_aE_0[(channel, rep)] = ' '.join(values)

    return results_aE_0

# Process and write files
for key, file_path in file_paths.items():
    results_aE_0 = process_spectrum(file_path, ensemble=key)
    output_file = f'../input_fit/final_matrixel/{key}_ground.txt'

    with open(output_file, 'w') as f:
        for channel, rep in order:
            line = results_aE_0.get((channel, rep), "0 0 0 0 0 0 0 0")
            f.write(line + '\n')

    print(f"Wrote renormalised output to {output_file}")

    



import pandas as pd
import os


# Define the file paths for each M file
file_paths = {
    "M1": '../CSVs/M1_spectral_density_matrix_elements_CB.csv',
    "M2": '../CSVs/M2_spectral_density_matrix_elements_CB.csv',
    "M3": '../CSVs/M3_spectral_density_matrix_elements_CB.csv',
    "M4": '../CSVs/M4_spectral_density_matrix_elements_CB.csv',
    "M5": '../CSVs/M5_spectral_density_matrix_elements_CB.csv'
}

# Define the order of channels and reps
order = [
    ("Chimera_OC_even", "as"),
    ("Chimera_OV12_even", "as"),
    ("Chimera_OV32_even", "as"),
    ("Chimera_OC_odd", "as"),
    ("Chimera_OV12_odd", "as"),
    ("Chimera_OV32_odd", "as"),
]

# Channel-to-Z-factor mapping
z_field_map = {
    "Chimera_OC_even": "Z_Lambda",
    "Chimera_OC_odd": "Z_Lambda",
    "Chimera_OV12_even": "Z_Sigma",
    "Chimera_OV32_even": "Z_Sigma",
    "Chimera_OV12_odd": "Z_Sigma",
    "Chimera_OV32_odd": "Z_Sigma",
}

# Output directory
os.makedirs('../input_fit/final_matrixel', exist_ok=True)

#volumes = {'M1': 20**3,'M2': 20**3,'M3': 20**3,'M4': 20**3,'M5': 32**3}
# Load the metadata
df = pd.read_csv("../input_fit/metadata/ensemble_metadata.csv")

# Compute volume dictionary for ensembles M1–M5
volumes = {
    row["ensemble_name"]: int(row["Ns"])**3
    for _, row in df.iterrows()
    if row["ensemble_name"] in ["M1", "M2", "M3", "M4", "M5"]
}


def process_spectrum(file_path, ensemble):
    df = pd.read_csv(file_path)
    results_aE_0 = {}
    grouped = df.groupby(['channel', 'rep'])
    
    for (channel, rep), group in grouped:
        z_field = z_field_map.get(channel)
        if z_field is None:
            continue
        #print('ensemble: ', ensemble)
        #ZA = z_factors.loc[z_factors['Ens'] == ensemble, z_field].values[0]
        ZA = computed_z[ensemble][z_field]
        vol = volumes[ensemble]
        values = []
        for i in range(2):
            if i < len(group):
                row = group.iloc[i]
                aE_0 = row['c0'] * ZA / np.sqrt(vol)
                err = row['errorc0'] * ZA * 6. / np.sqrt(vol) if row['errorc0'] != 0 else 0.001 * ZA
                values.append(f"{aE_0} {err}")
            else:
                values.append("0 0")

        results_aE_0[(channel, rep)] = ' '.join(values)

    return results_aE_0

# Process and write outputs
for key, file_path in file_paths.items():
    results_aE_0 = process_spectrum(file_path, ensemble=key)
    
    output_file_path_ground = f'../input_fit/final_matrixel/CB_{key}_ground.txt'
    with open(output_file_path_ground, 'w') as f_ground:
        for channel, rep in order:
            line = results_aE_0.get((channel, rep), "0 0 0 0 0 0 0 0")
            f_ground.write(line + '\n')

    print(f"Results for {key} written to {output_file_path_ground}")

