import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.legend_handler import HandlerTuple
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.cm as cm

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
    values_ground1.append(value * w0_M1)
    errors_ground1.append(error_diff * w0_M1)
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
    values_first1.append(value * w0_M1)
    errors_first1.append(error_diff * w0_M1)

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
    values_ground2.append(value * w0_M2)
    errors_ground2.append(error_diff * w0_M2)
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
    values_first2.append(value * w0_M2)
    errors_first2.append(error_diff * w0_M2)
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
    values_second2.append(value * w0_M2)
    errors_second2.append(error_diff * w0_M2)

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
    values_ground3.append(value * w0_M3)
    errors_ground3.append(error_diff * w0_M3)

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
    values_first3.append(value * w0_M3)
    errors_first3.append(error_diff * w0_M3)
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
    values_second3.append(value * w0_M3)
    errors_second3.append(error_diff * w0_M3)


# Read input file for M4_ground
with open('M4_ground.txt', 'r') as file:
    lines = file.readlines()

values_ground4 = []
errors_ground4 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_ground4.append(value * w0_M4)
    errors_ground4.append(error_diff * w0_M4)
print('\n')
# Read input file for M4_first
with open('M4_first.txt', 'r') as file:
    lines = file.readlines()

values_first4 = []
errors_first4 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_first4.append(value * w0_M4)
    errors_first4.append(error_diff * w0_M4)
print('\n')
# Read input file for M4_second
with open('M4_second.txt', 'r') as file:
    lines = file.readlines()

values_second4 = []
errors_second4 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_second4.append(value * w0_M4)
    errors_second4.append(error_diff * w0_M4)

print('\n\n')



# Read input file for M5_ground
with open('M5_ground.txt', 'r') as file:
    lines = file.readlines()

values_ground5 = []
errors_ground5 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_ground5.append(value * w0_M5)
    errors_ground5.append(error_diff * w0_M5)
print('\n')
# Read input file for M5_first
with open('M5_first.txt', 'r') as file:
    lines = file.readlines()

values_first5 = []
errors_first5 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_first5.append(value * w0_M5)
    errors_first5.append(error_diff * w0_M5)
print('\n')
# Read input file for M5_second
with open('M5_second.txt', 'r') as file:
    lines = file.readlines()

values_second5 = []
errors_second5 = []

# Process each line and store values and errors
for line in lines:
    value, stat, sys = process_line(line.strip())
    error_diff = math.sqrt(sys**2 + stat**2)
    values_second5.append(value * w0_M5)
    errors_second5.append(error_diff * w0_M5)

print('\n\n')
# Labels for x-axis
x_labels = ['PS', 'V', 'T', 'AV', 'AT', 'S', 'ps', 'v', 't', 'av', 'at', 's']

# Plotting
fig, ax = plt.subplots(figsize=(11, 6))


# Plot rectangles for each value and error_diff with a bit of spacing for M1_ground
spacing = 0.4
hatches = [r'//////', '++++']  # Define different hatch patterns

# Define the number of colors you want
num_colors = 6

# Generate evenly spaced values between 0 and 1
indices = np.linspace(0, 0.5, num_colors)

# Get the RGB colors from the Viridis colormap at these indices
viridis_colors = [cm.tab10(index) for index in indices]


for i in range(len(values_ground1)):
    if i < 6:
        color = viridis_colors[i]
    else:
        color = viridis_colors[i-6]
    #color = average_color1_4
    ax.add_patch(plt.Rectangle((i - 0.5 + spacing*i +0.2, values_ground1[i] - errors_ground1[i]), 0.4, 2*errors_ground1[i], color=color, alpha=0.4))
    ax.plot([i + spacing*i, i + spacing*i], [values_ground1[i] - errors_ground1[i], values_ground1[i] + errors_ground1[i]], color='red', alpha=0.0)
    
# Plot rectangles for each value and error_diff with a bit of spacing for M1_first
for i in range(len(values_first1)):
    if i < 6:
        color = viridis_colors[i]
    else:
        color = viridis_colors[i-6]
    ax.add_patch(plt.Rectangle((i - 0.5 + spacing*i +0.2, values_first1[i] - errors_first1[i]), 0.4, 2*errors_first1[i], color=color, alpha=0.4))
    ax.plot([i + spacing*i, i + spacing*i], [values_first1[i] - errors_first1[i], values_first1[i] + errors_first1[i]], color='orange', alpha=0.0)


# Plot rectangles for each value and error_diff with a bit of spacing for M2_ground
for i in range(len(values_ground2)):
    if i < 6:
        color = viridis_colors[i]
    else:
        color = viridis_colors[i-6]
    ax.add_patch(plt.Rectangle((i - 0.4 + spacing*i +0.2, values_ground2[i] - errors_ground2[i]), 0.4, 2*errors_ground2[i], color=color, alpha=0.6))
    ax.plot([i + spacing*i, i + spacing*i], [values_ground2[i] - errors_ground2[i], values_ground2[i] + errors_ground2[i]], color='blue', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M2_first
for i in range(len(values_first2)):
    if i < 6:
        color = viridis_colors[i]
    else:
        color = viridis_colors[i-6]
    ax.add_patch(plt.Rectangle((i - 0.4 + spacing*i +0.2, values_first2[i] - errors_first2[i]), 0.4, 2*errors_first2[i], color=color, alpha=0.6))
    ax.plot([i + spacing*i, i + spacing*i], [values_first2[i] - errors_first2[i], values_first2[i] + errors_first2[i]], color='green', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M2_second
for i in range(len(values_second2)):
    if i < 6:
        color = viridis_colors[i]
    else:
        color = viridis_colors[i-6]
    ax.add_patch(plt.Rectangle((i - 0.5 + spacing*i +0.2, values_second2[i] - errors_second2[i]), 0.4, 2*errors_second2[i], color=color, alpha=0.6))
    ax.plot([i + spacing*i, i + spacing*i], [values_second2[i] - errors_second2[i], values_second2[i] + errors_second2[i]], color='black', alpha=0.0)

#ax.add_patch(plt.Rectangle((- 0.5 , values_ground3[0] - errors_ground3[0]), 0.4, 2*errors_ground3[0], color=color3, alpha=0.7, label = '$N_t = 96$'))
# Plot rectangles for each value and error_diff with a bit of spacing for M3_ground
for i in range(len(values_ground3)):
    if i < 6:
        color = viridis_colors[i]
    else:
        color = viridis_colors[i-6]
    ax.add_patch(plt.Rectangle((i - 0.3 + spacing*i +0.2, values_ground3[i] - errors_ground3[i]), 0.4, 2*errors_ground3[i], color=color, alpha=0.85))
    ax.plot([i + spacing*i, i + spacing*i], [values_ground3[i] - errors_ground3[i], values_ground3[i] + errors_ground3[i]], color='blue', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M3_first
for i in range(len(values_first3)):
    if i < 6:
        color = viridis_colors[i]
    else:
        color = viridis_colors[i-6]
    ax.add_patch(plt.Rectangle((i - 0.3 + spacing*i +0.2, values_first3[i] - errors_first3[i]), 0.4, 2*errors_first3[i], color=color, alpha=0.85))
    ax.plot([i + spacing*i, i + spacing*i], [values_first3[i] - errors_first3[i], values_first3[i] + errors_first3[i]], color='green', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M3_second
for i in range(len(values_second3)):
    if i < 6:
        color = viridis_colors[i]
    else:
        color = viridis_colors[i-6]
    ax.add_patch(plt.Rectangle((i - 0.4 + spacing*i +0.2, values_second3[i] - errors_second3[i]), 0.4, 2*errors_second3[i], color=color, alpha=0.85))
    ax.plot([i + spacing*i, i + spacing*i], [values_second3[i] - errors_second3[i], values_second3[i] + errors_second3[i]], color='black', alpha=0.0)



######
import colorsys
# Function to adjust saturation of a color
def adjust_saturation(color, factor):
    # Convert color from RGB to HSL (Hue, Saturation, Lightness)
    h, l, s = colorsys.rgb_to_hls(*color[:3])

    # Increase saturation by multiplying with the factor
    s *= factor

    # Limit saturation to the range [0, 1]
    s = max(0.0, min(s, 1.0))

    # Convert back from HSL to RGB
    r, g, b = colorsys.hls_to_rgb(h, l, s)

    return (r, g, b, color[3])
# Plot rectangles for each value and error_diff with a bit of spacing for M4_ground
for i in range(len(values_ground3)):
    if i < 6:
        color = viridis_colors[i]
    else:
        color = viridis_colors[i-6]
    ax.add_patch(plt.Rectangle((i - 0.7 + spacing*i +0.2, values_ground4[i] - errors_ground4[i]), 0.15, 2*errors_ground4[i], color=adjust_saturation(color, 1.2), alpha=0.7, hatch=hatches[1],fill=False))
    ax.plot([i + spacing*i, i + spacing*i], [values_ground4[i] - errors_ground4[i], values_ground4[i] + errors_ground4[i]], color='blue', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M4_first
for i in range(len(values_first3)):
    if i < 6:
        color = viridis_colors[i]
    else:
        color = viridis_colors[i-6]
    ax.add_patch(plt.Rectangle((i - 0.7 + spacing*i +0.2, values_first4[i] - errors_first4[i]), 0.15, 2*errors_first4[i], color=adjust_saturation(color, 1.2), alpha=0.7, hatch=hatches[1],fill=False))
    ax.plot([i + spacing*i, i + spacing*i], [values_first4[i] - errors_first4[i], values_first4[i] + errors_first4[i]], color='green', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M4_second
for i in range(len(values_second3)):
    if i < 6:
        color = viridis_colors[i]
    else:
        color = viridis_colors[i-6]
    ax.add_patch(plt.Rectangle((i - 0.7 + spacing*i +0.2, values_second4[i] - errors_second4[i]), 0.15, 2*errors_second4[i], color=adjust_saturation(color, 1.2), alpha=0.7, hatch=hatches[1],fill=False))
    ax.plot([i + spacing*i, i + spacing*i], [values_second4[i] - errors_second4[i], values_second4[i] + errors_second4[i]], color='black', alpha=0.0)


######

# Plot rectangles for each value and error_diff with a bit of spacing for M4_ground
for i in range(len(values_ground3)):
    if i < 6:
        color = viridis_colors[i]
    else:
        color = viridis_colors[i-6]
    ax.add_patch(plt.Rectangle((i + 0.15 + spacing*i+0.2, values_ground5[i] - errors_ground5[i]), 0.15, 2*errors_ground5[i], color=adjust_saturation(color, 1.2), alpha=0.7, hatch=hatches[0],fill=False))
    ax.plot([i + spacing*i, i + spacing*i], [values_ground5[i] - errors_ground5[i], values_ground5[i] + errors_ground5[i]], color='blue', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M4_first
for i in range(len(values_first3)):
    if i < 6:
        color = viridis_colors[i]
    else:
        color = viridis_colors[i-6]
    ax.add_patch(plt.Rectangle((i + 0.15 + spacing*i+0.2, values_first5[i] - errors_first5[i]), 0.15, 2*errors_first5[i], color=adjust_saturation(color, 1.2), alpha=0.7, hatch=hatches[0],fill=False))
    ax.plot([i + spacing*i, i + spacing*i], [values_first5[i] - errors_first5[i], values_first5[i] + errors_first5[i]], color='green', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M4_second
for i in range(len(values_second3)):
    if i < 6:
        color = viridis_colors[i]
    else:
        color = viridis_colors[i-6]
    ax.add_patch(plt.Rectangle((i + 0.15 + spacing*i+0.2, values_second5[i] - errors_second5[i]), 0.15, 2*errors_second5[i], color=adjust_saturation(color, 1.2), alpha=0.7, hatch=hatches[0],fill=False))
    ax.plot([i + spacing*i, i + spacing*i], [values_second5[i] - errors_second5[i], values_second5[i] + errors_second5[i]], color='black', alpha=0.0)

# Set labels and title
#ax.set_xlabel('Channels', y=-0.80, fontsize=14)
ax.set_ylabel('$\hat{m}$', fontsize=16)
#ax.set_title('Mesonic spectrum, $am^{\\rm f}_0 = -0.71$, $am^{\\rm as}_0 = -1.01$', fontsize=16)
ax.set_title('Meson spectrum', fontsize=16)

# Set x-axis ticks and labels
ax.set_xticks([i + spacing*i for i in range(len(values_ground1))])
ax.set_xticklabels(x_labels, ha='center')



# Create custom legend handles with color patches for each group
legend_handles = [
    mpatches.Patch(color='black', alpha=0.4),
    mpatches.Patch(color='black', alpha=0.6),
    mpatches.Patch(color='black', alpha=0.85),
    mpatches.Patch(color='black', alpha=0.7,fill=False),
    mpatches.Patch(color='black', alpha=0.7,fill=False),
]

for i, hatch in enumerate(hatches):
    legend_handles[-(i+1)].set_hatch(hatch)

# Plot the legend with custom handles and labels
#plt.legend(handles=legend_handles, labels=legend_labels, handlelength=1.5, handletextpad=3.4)
plt.legend([(legend_handles[0]), (legend_handles[1]), (legend_handles[2]), (legend_handles[3]), (legend_handles[4])], ['M1', 'M2', 'M3', 'M4', 'M5'],
               handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=3.0)




plt.ylim(0.79,3.70)

plt.tight_layout()
plt.savefig(output_file, format='pdf')
#plt.show()

