import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# Activating text rendering by LaTeX
plt.style.use("paperdraft.mplstyle")

# Define the output file name
output_file = '../../plots/final_spectrum_M5.pdf'


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

    if min_error_index < len(errors) - 1:
        error_diff = math.sqrt(errors[min_error_index]**2 + errors[min_error_index+1]**2)
    else:
        error_diff = errors[min_error_index]
    
    stat_err = errors[min_error_index]
    return stat_measurement, stat_err, max_diff

# Read input file for M5_ground
with open('M5_ground.txt', 'r') as file:
    lines = file.readlines()

values_ground4 = []
errors_ground4 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_ground4.append(value)
    errors_ground4.append(error_diff)

# Read input file for M5_first
with open('M5_first.txt', 'r') as file:
    lines = file.readlines()

values_first4 = []
errors_first4 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_first4.append(value)
    errors_first4.append(error_diff)

# Read input file for M5_second
with open('M5_second.txt', 'r') as file:
    lines = file.readlines()

values_second4 = []
errors_second4 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_second4.append(value)
    errors_second4.append(error_diff)

# Labels for x-axis
x_labels = ['PS', 'V', 'T', 'AV', 'AT', 'S', 'ps', 'v', 't', 'av', 'at', 's']

# Plotting
fig, ax = plt.subplots(figsize=(10, 6))

color1  = cm.Oranges((20 + 12 + 6) / (12 + 12 + 12 + 12 + 12 + 12 + 12) )
color2  = cm.Oranges((20 + 12 + 12 + 12) / (12 + 12 + 12 + 12 + 12 + 12 + 12) )
color3  = cm.Oranges((20 + 12 + 12 + 12 + 12 + 12) / (12 + 12 + 12 + 12 + 12 + 12 + 12) )

# Plot rectangles for each value and error_diff with a bit of spacing for M5_ground
spacing = 0.2
ax.add_patch(plt.Rectangle(( - 0.5 , values_ground4[0] - errors_ground4[0]), 1, 2*errors_ground4[0], color=color1, alpha=0.7, label='$E_0$'))
for i in range(1,len(values_ground4)):
    ax.add_patch(plt.Rectangle((i - 0.5 + spacing*i, values_ground4[i] - errors_ground4[i]), 1, 2*errors_ground4[i], color=color1, alpha=0.7))
    ax.plot([i + spacing*i, i + spacing*i], [values_ground4[i] - errors_ground4[i], values_ground4[i] + errors_ground4[i]], color='red', alpha=0.0)

ax.add_patch(plt.Rectangle(( - 0.5 , values_first4[0] - errors_first4[0]), 1, 2*errors_first4[0], color=color2, alpha=0.7, label='$E_1$'))
# Plot rectangles for each value and error_diff with a bit of spacing for M5_first
for i in range(1,len(values_first4)):
    ax.add_patch(plt.Rectangle((i - 0.5 + spacing*i, values_first4[i] - errors_first4[i]), 1, 2*errors_first4[i], color=color2, alpha=0.7))
    ax.plot([i + spacing*i, i + spacing*i], [values_first4[i] - errors_first4[i], values_first4[i] + errors_first4[i]], color='orange', alpha=0.0)

ax.add_patch(plt.Rectangle(( - 0.5 , values_second4[0] - errors_second4[0]), 1, 2*errors_second4[0], color=color3, alpha=0.7, label='$E_2$'))
# Plot rectangles for each value and error_diff with a bit of spacing for M5_second
for i in range(1,len(values_second4)):
    ax.add_patch(plt.Rectangle((i - 0.5 + spacing*i, values_second4[i] - errors_second4[i]), 1, 2*errors_second4[i], color=color3, alpha=0.7))
    ax.plot([i + spacing*i, i + spacing*i], [values_second4[i] - errors_second4[i], values_second4[i] + errors_second4[i]], color='black', alpha=0.0)

# Set labels and title
#ax.set_xlabel('Channels', y=-0.80, fontsize=14)
ax.set_ylabel('Energy levels $[a^{-1}]$', fontsize=14)
ax.set_title('Mesonic spectrum, $am^{\\rm f}_0 = -0.72$, $am^{\\rm as}_0 = -1.01$', fontsize=16)

# Set x-axis ticks and labels
ax.set_xticks([i + spacing*i for i in range(len(values_ground4))])
ax.set_xticklabels(x_labels, ha='center')

# Add legend
ax.legend()

plt.ylim(0.29, 1.20)

plt.tight_layout()
plt.savefig(output_file, format='pdf')
#plt.show()

