import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.legend_handler import HandlerTuple
import matplotlib.patches as mpatches



# Activating text rendering by LaTeX
plt.style.use("paperdraft.mplstyle")

# Define the output file name
output_file = '../../plots/final_spectrum_M1_M2_M3.pdf'


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
    print('Meas: ', stat_measurement, '\tstat: ' , stat_err, '\tsys: ', max_diff, '\tQuad: ', math.sqrt(max_diff**2 + stat_err**2))
    return stat_measurement, stat_err, max_diff

# Read input file for M1_ground
with open('M1_ground.txt', 'r') as file:
    lines = file.readlines()

values_ground1 = []
errors_ground1 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_ground1.append(value)
    errors_ground1.append(error_diff)
print('\n')

# Read input file for M1_first
with open('M1_first.txt', 'r') as file:
    lines = file.readlines()

values_first1 = []
errors_first1 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_first1.append(value)
    errors_first1.append(error_diff)

print('\n\n')

# Read input file for M2_ground
with open('M2_ground.txt', 'r') as file:
    lines = file.readlines()

values_ground2 = []
errors_ground2 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_ground2.append(value)
    errors_ground2.append(error_diff)
print('\n')
# Read input file for M2_first
with open('M2_first.txt', 'r') as file:
    lines = file.readlines()

values_first2 = []
errors_first2 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_first2.append(value)
    errors_first2.append(error_diff)
print('\n')
# Read input file for M2_second
with open('M2_second.txt', 'r') as file:
    lines = file.readlines()

values_second2 = []
errors_second2 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_second2.append(value)
    errors_second2.append(error_diff)

print('\n\n')

# Read input file for M3_ground
with open('M3_ground.txt', 'r') as file:
    lines = file.readlines()

values_ground3 = []
errors_ground3 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_ground3.append(value)
    errors_ground3.append(error_diff)

print('\n')
# Read input file for M3_first
with open('M3_first.txt', 'r') as file:
    lines = file.readlines()

values_first3 = []
errors_first3 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_first3.append(value)
    errors_first3.append(error_diff)
print('\n')
# Read input file for M3_second
with open('M3_second.txt', 'r') as file:
    lines = file.readlines()

values_second3 = []
errors_second3 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_second3.append(value)
    errors_second3.append(error_diff)

# Labels for x-axis
x_labels = ['PS', 'V', 'T', 'AV', 'AT', 'S', 'ps', 'v', 't', 'av', 'at', 's']

# Plotting
fig, ax = plt.subplots(figsize=(10, 6))

color1  = cm.Oranges((20 +len(values_ground1)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )
color2  = cm.Oranges(( 20 +len(values_ground1) + len(values_first1)+ len(values_ground2)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )
color3 = cm.Oranges((20 +len(values_ground1) + len(values_first1)+ len(values_ground2)+len(values_first2)+len(values_second2)+len(values_first3)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )
color4  = cm.Blues((20 +len(values_ground1)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )
color5  = cm.Blues(( 20 +len(values_ground1) + len(values_first1)+ len(values_ground2)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )
color6 = cm.Blues((20 +len(values_ground1) + len(values_first1)+ len(values_ground2)+len(values_first2)+len(values_second2)+len(values_first3)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )


# Plot rectangles for each value and error_diff with a bit of spacing for M1_ground
spacing = 0.4
color = color1
ax.add_patch(plt.Rectangle((- 0.5 , values_ground1[0] - errors_ground1[0]), 1, 2*errors_ground1[0], color=color, alpha=0.7, label = '$N_t = 48$'))
for i in range(len(values_ground1)):
    if i < 6:
        color = color1
    else:
        color = color4
    ax.add_patch(plt.Rectangle((i - 0.5 + spacing*i, values_ground1[i] - errors_ground1[i]), 1, 2*errors_ground1[i], color=color, alpha=0.7))
    ax.plot([i + spacing*i, i + spacing*i], [values_ground1[i] - errors_ground1[i], values_ground1[i] + errors_ground1[i]], color='red', alpha=0.0)
    
# Plot rectangles for each value and error_diff with a bit of spacing for M1_first
for i in range(len(values_first1)):
    if i < 6:
        color = color1
    else:
        color = color4
    ax.add_patch(plt.Rectangle((i - 0.5 + spacing*i, values_first1[i] - errors_first1[i]), 1, 2*errors_first1[i], color=color, alpha=0.7))
    ax.plot([i + spacing*i, i + spacing*i], [values_first1[i] - errors_first1[i], values_first1[i] + errors_first1[i]], color='orange', alpha=0.0)

ax.add_patch(plt.Rectangle((- 0.5 , values_ground2[0] - errors_ground2[0]), 1, 2*errors_ground2[0], color=color2, alpha=0.7, label = '$N_t = 64$'))
# Plot rectangles for each value and error_diff with a bit of spacing for M2_ground
for i in range(len(values_ground2)):
    if i < 6:
        color = color2
    else:
        color = color5
    ax.add_patch(plt.Rectangle((i - 0.4 + spacing*i, values_ground2[i] - errors_ground2[i]), 1, 2*errors_ground2[i], color=color, alpha=0.7))
    ax.plot([i + spacing*i, i + spacing*i], [values_ground2[i] - errors_ground2[i], values_ground2[i] + errors_ground2[i]], color='blue', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M2_first
for i in range(len(values_first2)):
    if i < 6:
        color = color2
    else:
        color = color5
    ax.add_patch(plt.Rectangle((i - 0.4 + spacing*i, values_first2[i] - errors_first2[i]), 1, 2*errors_first2[i], color=color, alpha=0.7))
    ax.plot([i + spacing*i, i + spacing*i], [values_first2[i] - errors_first2[i], values_first2[i] + errors_first2[i]], color='green', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M2_second
for i in range(len(values_second2)):
    if i < 6:
        color = color2
    else:
        color = color5
    ax.add_patch(plt.Rectangle((i - 0.5 + spacing*i, values_second2[i] - errors_second2[i]), 1, 2*errors_second2[i], color=color, alpha=0.7))
    ax.plot([i + spacing*i, i + spacing*i], [values_second2[i] - errors_second2[i], values_second2[i] + errors_second2[i]], color='black', alpha=0.0)

ax.add_patch(plt.Rectangle((- 0.5 , values_ground3[0] - errors_ground3[0]), 1, 2*errors_ground3[0], color=color3, alpha=0.7, label = '$N_t = 96$'))
# Plot rectangles for each value and error_diff with a bit of spacing for M3_ground
for i in range(len(values_ground3)):
    if i < 6:
        color = color3
    else:
        color = color6
    ax.add_patch(plt.Rectangle((i - 0.3 + spacing*i, values_ground3[i] - errors_ground3[i]), 1, 2*errors_ground3[i], color=color, alpha=0.7))
    ax.plot([i + spacing*i, i + spacing*i], [values_ground3[i] - errors_ground3[i], values_ground3[i] + errors_ground3[i]], color='blue', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M3_first
for i in range(len(values_first3)):
    if i < 6:
        color = color3
    else:
        color = color6
    ax.add_patch(plt.Rectangle((i - 0.3 + spacing*i, values_first3[i] - errors_first3[i]), 1, 2*errors_first3[i], color=color, alpha=0.7))
    ax.plot([i + spacing*i, i + spacing*i], [values_first3[i] - errors_first3[i], values_first3[i] + errors_first3[i]], color='green', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M3_second
for i in range(len(values_second3)):
    if i < 6:
        color = color3
    else:
        color = color6
    ax.add_patch(plt.Rectangle((i - 0.4 + spacing*i, values_second3[i] - errors_second3[i]), 1, 2*errors_second3[i], color=color, alpha=0.7))
    ax.plot([i + spacing*i, i + spacing*i], [values_second3[i] - errors_second3[i], values_second3[i] + errors_second3[i]], color='black', alpha=0.0)

# Set labels and title
#ax.set_xlabel('Channels', y=-0.80, fontsize=14)
ax.set_ylabel('Energy levels $[a^{-1}]$', fontsize=14)
ax.set_title('Mesonic spectrum, $am^{\\rm f}_0 = -0.71$, $am^{\\rm as}_0 = -1.01$', fontsize=16)

# Set x-axis ticks and labels
ax.set_xticks([i + spacing*i for i in range(len(values_ground1))])
ax.set_xticklabels(x_labels, ha='center')

# Create custom legend handles with color patches for each group
legend_handles = [
    mpatches.Patch(color=color1, label='\t', alpha=0.7),
    mpatches.Patch(color=color4, alpha=0.7),
    mpatches.Patch(color=color2, label='\t', alpha=0.7),
    mpatches.Patch(color=color5, alpha=0.7),
    mpatches.Patch(color=color3, label='\t ', alpha=0.7),
    mpatches.Patch(color=color6, alpha=0.7)
]

'''
# Add legend
ax.legend(handles=legend_handles, handlelength=1.5, handletextpad=3.4)

# Insert text at position (x, y)
plt.text(x=-0.20, y=1.39, s='$N_t = 48$', fontsize=14, color='black', alpha=1.0)
# Insert text at position (x, y)
plt.text(x=-0.20, y=1.255, s='$N_t = 64$', fontsize=14, color='black')
# Insert text at position (x, y)
plt.text(x=-0.20, y=1.115, s='$N_t = 96$', fontsize=14, color='black')

# Iterate over legend handles and set text color to black
for text in ax.get_legend().get_texts():
    text.set_color('black')
'''
plt.ylim(0.33,1.50)

plt.tight_layout()
plt.savefig(output_file, format='pdf')
#plt.show()

