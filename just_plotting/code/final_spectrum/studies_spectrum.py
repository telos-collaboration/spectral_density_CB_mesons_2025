import os
import math
import matplotlib.pyplot as plt
import numpy as np

w0_M1 = 2.5100
w0_M2 = 2.5300
w0_M3 = 2.5170
w0_M4 = 2.3557
w0_M5 = 2.6927

# Activating text rendering by LaTeX
plt.style.use("paperdraft.mplstyle")

# Define the output file name
output_file = '../../../plots/final_spectrum_MN3.pdf'


def process_line(line):
    # Split the line into measurements and errors
    data = [float(x) for x in line.split()]
    
    # Extract measurements and errors
    measurements = data[::2]
    errors = data[1::2]

    # Find the index of the minimum error
    min_error_index = errors.index(min(errors))

    # Calculate the stat measurement (smallest error measurement)
    stat_measurement = measurements[min_error_index]

    # Calculate sys (maximum difference between measurements)
    diffs = [abs(measurements[i] - measurements[j]) for i in range(len(measurements)) for j in range(i+1, len(measurements))]

    max_diff = max(diffs)
    
    stat_err = errors[min_error_index]
    #print('Meas: ', stat_measurement, '\tstat: ' , stat_err, '\tsys: ', max_diff, '\tQuad: ', math.sqrt(max_diff**2 + stat_err**2))
    return stat_measurement, stat_err, max_diff

# Function to read input file and process data
def read_data(filename):
    if not os.path.exists(filename):
        return None, None
    
    with open(filename, 'r') as file:
        lines = file.readlines()

    values = []
    errors = []

    # Process each line and store values and errors
    for line in lines:
        value, stat, sys = process_line(line.strip())
        error_diff = math.sqrt(sys**2 + stat**2)
        values.append(value)
        errors.append(error_diff)
    
    return values, errors

# Read input data for ground state
values_ground = []
errors_ground = []
for i in range(1, 6):
    filename = f'M{i}_ground.txt'
    values, errors = read_data(filename)
    if values is not None and errors is not None:
        values_ground.append(values)
        errors_ground.append(errors)

# Read input data for first excited state
values_first = []
errors_first = []
for i in range(1, 6):
    filename = f'M{i}_first.txt'
    values, errors = read_data(filename)
    if values is not None and errors is not None:
        values_first.append(values)
        errors_first.append(errors)

# Read input data for second excited state
values_second = []
errors_second = []
for i in range(1, 6):
    filename = f'M{i}_second.txt'
    values, errors = read_data(filename)
    if values is not None and errors is not None:
        values_second.append(values)
        errors_second.append(errors)

# Define colors for each dataset
colors = ['blue', 'green', 'red', 'purple', 'orange']

# Plot all the data together
plt.figure(figsize=(8, 6))

# Plot ground state data
for i in range(len(values_ground)):
    plt.errorbar(values_ground[i][0], 
                 values_ground[i][6], 
                 xerr=errors_ground[i][0], 
                 yerr=errors_ground[i][6], 
                 fmt='o', 
                 capsize=5, 
                 label=f'M{i+1} (Ground)', 
                 color=colors[i], 
                 elinewidth=1.4, 
                 markersize=7.0)

# Plot first excited state data
for i in range(len(values_first)):
    plt.errorbar(values_first[i][0], 
                 values_first[i][6], 
                 xerr=errors_first[i][0], 
                 yerr=errors_first[i][6], 
                 fmt='s', 
                 capsize=5, 
                 label=f'M{i+1} (First)', 
                 color=colors[i], 
                 elinewidth=1.4, 
                 markersize=7.0)

# Plot second excited state data if available
for i in range(1,2):
    plt.errorbar(values_second[i][0], 
                 values_second[i][6], 
                 xerr=errors_second[i][0], 
                 yerr=errors_second[i][6], 
                 fmt='^', 
                 capsize=5, 
                 label=f'M{i+1} (Second)', 
                 color=colors[i], 
                 elinewidth=1.4, 
                 markersize=7.0)

plt.xlabel('$m_{\mathrm PS}$', fontsize=14)
plt.ylabel('$m_{\mathrm ps}$', fontsize=14)
#plt.title('Comparison of Ground, First Excited, and Second Excited States', fontsize=14)
plt.grid(True, linestyle='--')
plt.legend(loc='best', fontsize=12)
#plt.xlim(0.83, 0.974, 0.2)

for i in range(len(values_first)):
    print(f'M{i+1} Ratios ground: ', values_ground[i][0]/values_ground[i][6], '+-', np.sqrt((errors_ground[i][0]/values_ground[i][6])**2 + (values_ground[i][0]*errors_ground[i][6]/values_ground[i][6]**2)**2))
print('\n')
for i in range(len(values_first)):
    print(f'M{i+1} Ratios first: ', values_first[i][0]/values_first[i][6], '+-', np.sqrt((errors_first[i][0]/values_first[i][6])**2 + (values_first[i][0]*errors_first[i][6]/values_first[i][6]**2)**2))
print('\n')
for i in range(1,2):
    print(f'M{i+1} Ratios second: ', values_second[i][0]/values_second[i][6], '+-', np.sqrt((errors_second[i][0]/values_second[i][6])**2 + (values_second[i][0]*errors_second[i][6]/values_second[i][6]**2)**2))


# Save the plot to a PDF
plt.savefig('../../../plots/PS_ground_first_second_combined.pdf')

# Show the plot
#plt.show()
