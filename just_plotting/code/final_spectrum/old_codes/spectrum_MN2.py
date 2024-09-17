import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.legend_handler import HandlerTuple
import matplotlib.patches as mpatches

w0_M1 = 2.5100
w0_M2 = 2.5300
w0_M3 = 2.5170
w0_M4 = 2.3557
w0_M5 = 2.6927

# Activating text rendering by LaTeX
plt.style.use("paperdraft.mplstyle")

# Define the output file name
output_file = '../../plots/final_spectrum_MN2.pdf'


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

color1  = cm.Oranges((20 +len(values_ground1)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )
color2  = cm.Oranges(( 20 +len(values_ground1) + len(values_first1)+ len(values_ground2)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )
color3 = cm.Oranges((20 +len(values_ground1) + len(values_first1)+ len(values_ground2)+len(values_first2)+len(values_second2)+len(values_first3)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )
color7 = cm.cividis(( len(values_ground1) + len(values_first1 ) + 10 + len(values_ground2)+len(values_first2)+len(values_second2)+len(values_first3)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )
color8 = cm.cividis((len(values_ground1) + len(values_first1) - 5 + len(values_ground2)+len(values_first2)+len(values_second2)+len(values_first3)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )

color4  = cm.plasma((20 +len(values_ground1)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )
color5  = cm.plasma(( 20 +len(values_ground1) + len(values_first1)+ len(values_ground2)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )
color6 = cm.plasma((20 +len(values_ground1) + len(values_first1)+ len(values_ground2)+len(values_first2)+len(values_second2)+len(values_first3)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )
color9 = cm.magma(( len(values_ground1) + len(values_first1 ) - 30 + len(values_ground2)+len(values_first2)+len(values_second2)+len(values_first3)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )
color10 = cm.magma((len(values_ground1) + len(values_first1) + 5 + len(values_ground2)+len(values_first2)+len(values_second2)+len(values_first3)) / (len(values_ground1) + len(values_first1) + len(values_ground2) + len(values_first2)+len(values_ground3) + len(values_first3)+len(values_second3)) )


import matplotlib.colors as mcolors

def average_color(color1, color2):
    rgba1 = mcolors.to_rgba(color1)
    rgba2 = mcolors.to_rgba(color2)
    avg_rgba = [0.7*(c1 + c2) / 3 for c1, c2 in zip(rgba1, rgba2)]
    return mcolors.to_hex(avg_rgba)
    

def average_color2(color1, color2):
    rgba1 = mcolors.to_rgba(color1)
    rgba2 = mcolors.to_rgba(color2)
    avg_rgba = [0.8*(c1 + c2) / 3 for c1, c2 in zip(rgba1, rgba2)]
    return mcolors.to_hex(avg_rgba)    

def average_color3(color1, color2):
    rgba1 = mcolors.to_rgba(color1)
    rgba2 = mcolors.to_rgba(color2)
    avg_rgba = [1.0*(c1 + c2) / 3 for c1, c2 in zip(rgba1, rgba2)]
    return mcolors.to_hex(avg_rgba)  

def average_color4(color1, color2):
    rgba1 = mcolors.to_rgba(color1)
    rgba2 = mcolors.to_rgba(color2)
    avg_rgba = [7*(c1 + c2) / 18 for c1, c2 in zip(rgba1, rgba2)]
    return mcolors.to_hex(avg_rgba)  

def average_color5(color1, color2):
    rgba1 = mcolors.to_rgba(color1)
    rgba2 = mcolors.to_rgba(color2)
    avg_rgba = [8*(c1 + c2) / 18 for c1, c2 in zip(rgba1, rgba2)]
    return mcolors.to_hex(avg_rgba)   
    
def average_color6(color1, color2):
    rgba1 = mcolors.to_rgba(color1)
    rgba2 = mcolors.to_rgba(color2)
    avg_rgba = [9*(c1 + c2) / 18  for c1, c2 in zip(rgba1, rgba2)]
    return mcolors.to_hex(avg_rgba)  

# Calculate average color
average_color1_4 = average_color(color1, color4)
average_color1_4_2 = average_color2(color1, color4)
average_color1_4_3 = average_color3(color1, color4)
average_color1_4_4 = average_color4(color1, color4)
average_color1_4_5 = average_color5(color1, color4)
average_color1_4_6 = average_color6(color1, color4)
avg_1_4 = [average_color1_4, average_color1_4_2, average_color1_4_3, average_color1_4_4, average_color1_4_5, average_color1_4_6]


average_color2_5 = average_color(color2, color5)
average_color2_5_2 = average_color2(color2, color5)
average_color2_5_3 = average_color3(color2, color5)
average_color2_5_4 = average_color4(color2, color5)
average_color2_5_5 = average_color5(color2, color5)
average_color2_5_6 = average_color6(color2, color5)
avg_2_5 = [average_color2_5, average_color2_5_2, average_color2_5_3, average_color2_5_4, average_color2_5_5, average_color2_5_6]

average_color3_6 = average_color(color3, color6)
average_color3_6_2 = average_color2(color3, color6)
average_color3_6_3 = average_color3(color3, color6)
average_color3_6_4 = average_color4(color3, color6)
average_color3_6_5 = average_color5(color3, color6)
average_color3_6_6 = average_color6(color3, color6)
avg_3_6 = [average_color3_6, average_color3_6_2, average_color3_6_3, average_color3_6_4, average_color3_6_5, average_color3_6_6]



average_color7_9 = average_color(color7, color9)
average_color7_9_2 = average_color2(color7, color9)
average_color7_9_3 = average_color3(color7, color9)
average_color7_9_4 = average_color4(color7, color9)
average_color7_9_5 = average_color5(color7, color9)
average_color7_9_6 = average_color6(color7, color9)
avg_7_9 = [average_color7_9, average_color7_9_2, average_color7_9_3, average_color7_9_4, average_color7_9_5, average_color7_9_6]

average_color8_10 = average_color(color8, color10)
average_color8_10_2 = average_color2(color8, color10)
average_color8_10_3 = average_color3(color8, color10)
average_color8_10_4 = average_color4(color8, color10)
average_color8_10_5 = average_color5(color8, color10)
average_color8_10_6 = average_color6(color8, color10)
avg_8_10 = [average_color8_10, average_color8_10_2, average_color8_10_3, average_color8_10_4, average_color8_10_5, average_color8_10_6]

# Plot rectangles for each value and error_diff with a bit of spacing for M1_ground
spacing = 0.4
hatches = [r'/////', 'O', '..', '--', '++', 'xx']  # Define different hatch patterns
color = average_color1_4
ax.add_patch(plt.Rectangle((- 0.5 +0.2 , values_ground1[0] - errors_ground1[0]), 0.4, 2*errors_ground1[0], color=color, alpha=0.7, label = '$N_t = 48$'))
for i in range(len(values_ground1)):
    if i < 6:
        color = avg_1_4[i]
    else:
        color = avg_1_4[i-6]
    #color = average_color1_4
    ax.add_patch(plt.Rectangle((i - 0.5 + spacing*i +0.2, values_ground1[i] - errors_ground1[i]), 0.4, 2*errors_ground1[i], color=color, alpha=0.7, hatch=hatches[i%len(hatches)]))
    ax.plot([i + spacing*i, i + spacing*i], [values_ground1[i] - errors_ground1[i], values_ground1[i] + errors_ground1[i]], color='red', alpha=0.0)
    
# Plot rectangles for each value and error_diff with a bit of spacing for M1_first
for i in range(len(values_first1)):
    if i < 6:
        color = avg_1_4[i]
    else:
        color = avg_1_4[i-6]
    ax.add_patch(plt.Rectangle((i - 0.5 + spacing*i +0.2, values_first1[i] - errors_first1[i]), 0.4, 2*errors_first1[i], color=color, alpha=0.7, hatch=hatches[i%len(hatches)]))
    ax.plot([i + spacing*i, i + spacing*i], [values_first1[i] - errors_first1[i], values_first1[i] + errors_first1[i]], color='orange', alpha=0.0)

ax.add_patch(plt.Rectangle((- 0.5 +0.2 , values_ground2[0] - errors_ground2[0]), 0.4, 2*errors_ground2[0], color=color2, alpha=0.7, label = '$N_t = 64$'))
# Plot rectangles for each value and error_diff with a bit of spacing for M2_ground
for i in range(len(values_ground2)):
    if i < 6:
        color = avg_2_5[i]
    else:
        color = avg_2_5[i-6]
    ax.add_patch(plt.Rectangle((i - 0.4 + spacing*i +0.2, values_ground2[i] - errors_ground2[i]), 0.4, 2*errors_ground2[i], color=color, alpha=0.7, hatch=hatches[i%len(hatches)]))
    ax.plot([i + spacing*i, i + spacing*i], [values_ground2[i] - errors_ground2[i], values_ground2[i] + errors_ground2[i]], color='blue', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M2_first
for i in range(len(values_first2)):
    if i < 6:
        color = avg_2_5[i]
    else:
        color = avg_2_5[i-6]
    ax.add_patch(plt.Rectangle((i - 0.4 + spacing*i +0.2, values_first2[i] - errors_first2[i]), 0.4, 2*errors_first2[i], color=color, alpha=0.7, hatch=hatches[i%len(hatches)]))
    ax.plot([i + spacing*i, i + spacing*i], [values_first2[i] - errors_first2[i], values_first2[i] + errors_first2[i]], color='green', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M2_second
for i in range(len(values_second2)):
    if i < 6:
        color = avg_2_5[i]
    else:
        color = avg_2_5[i-6]
    ax.add_patch(plt.Rectangle((i - 0.5 + spacing*i +0.2, values_second2[i] - errors_second2[i]), 0.4, 2*errors_second2[i], color=color, alpha=0.7, hatch=hatches[i%len(hatches)]))
    ax.plot([i + spacing*i, i + spacing*i], [values_second2[i] - errors_second2[i], values_second2[i] + errors_second2[i]], color='black', alpha=0.0)

#ax.add_patch(plt.Rectangle((- 0.5 , values_ground3[0] - errors_ground3[0]), 0.4, 2*errors_ground3[0], color=color3, alpha=0.7, label = '$N_t = 96$'))
# Plot rectangles for each value and error_diff with a bit of spacing for M3_ground
for i in range(len(values_ground3)):
    if i < 6:
        color = avg_3_6[i]
    else:
        color = avg_3_6[i-6]
    ax.add_patch(plt.Rectangle((i - 0.3 + spacing*i +0.2, values_ground3[i] - errors_ground3[i]), 0.4, 2*errors_ground3[i], color=color, alpha=0.7, hatch=hatches[i%len(hatches)]))
    ax.plot([i + spacing*i, i + spacing*i], [values_ground3[i] - errors_ground3[i], values_ground3[i] + errors_ground3[i]], color='blue', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M3_first
for i in range(len(values_first3)):
    if i < 6:
        color = avg_3_6[i]
    else:
        color = avg_3_6[i-6]
    ax.add_patch(plt.Rectangle((i - 0.3 + spacing*i +0.2, values_first3[i] - errors_first3[i]), 0.4, 2*errors_first3[i], color=color, alpha=0.7, hatch=hatches[i%len(hatches)]))
    ax.plot([i + spacing*i, i + spacing*i], [values_first3[i] - errors_first3[i], values_first3[i] + errors_first3[i]], color='green', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M3_second
for i in range(len(values_second3)):
    if i < 6:
        color = avg_3_6[i]
    else:
        color = avg_3_6[i-6]
    ax.add_patch(plt.Rectangle((i - 0.4 + spacing*i +0.2, values_second3[i] - errors_second3[i]), 0.4, 2*errors_second3[i], color=color, alpha=0.7, hatch=hatches[i%len(hatches)]))
    ax.plot([i + spacing*i, i + spacing*i], [values_second3[i] - errors_second3[i], values_second3[i] + errors_second3[i]], color='black', alpha=0.0)



######


ax.add_patch(plt.Rectangle((- 0.7 +0.2 , values_ground4[0] - errors_ground4[0]), 0.15, 2*errors_ground4[0], color=color7, alpha=0.7, label = '$N_t = 96$'))
# Plot rectangles for each value and error_diff with a bit of spacing for M4_ground
for i in range(len(values_ground3)):
    if i < 6:
        color = avg_7_9[i]
    else:
        color = avg_7_9[i-6]
    ax.add_patch(plt.Rectangle((i - 0.7 + spacing*i +0.2, values_ground4[i] - errors_ground4[i]), 0.15, 2*errors_ground4[i], color=color, alpha=0.7, hatch=hatches[i%len(hatches)]))
    ax.plot([i + spacing*i, i + spacing*i], [values_ground4[i] - errors_ground4[i], values_ground4[i] + errors_ground4[i]], color='blue', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M4_first
for i in range(len(values_first3)):
    if i < 6:
        color = avg_7_9[i]
    else:
        color = avg_7_9[i-6]
    ax.add_patch(plt.Rectangle((i - 0.7 + spacing*i +0.2, values_first4[i] - errors_first4[i]), 0.15, 2*errors_first4[i], color=color, alpha=0.7, hatch=hatches[i%len(hatches)]))
    ax.plot([i + spacing*i, i + spacing*i], [values_first4[i] - errors_first4[i], values_first4[i] + errors_first4[i]], color='green', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M4_second
for i in range(len(values_second3)):
    if i < 6:
        color = avg_7_9[i]
    else:
        color = avg_7_9[i-6]
    ax.add_patch(plt.Rectangle((i - 0.7 + spacing*i +0.2, values_second4[i] - errors_second4[i]), 0.15, 2*errors_second4[i], color=color, alpha=0.7, hatch=hatches[i%len(hatches)]))
    ax.plot([i + spacing*i, i + spacing*i], [values_second4[i] - errors_second4[i], values_second4[i] + errors_second4[i]], color='black', alpha=0.0)


######

ax.add_patch(plt.Rectangle((+ 0.15 +0.2, values_ground5[0] - errors_ground5[0]), 0.15, 2*errors_ground5[0], color=color8, alpha=0.7, label = '$N_t = 96$'))
# Plot rectangles for each value and error_diff with a bit of spacing for M4_ground
for i in range(len(values_ground3)):
    if i < 6:
        color = avg_8_10[i]
    else:
        color = avg_8_10[i-6]
    ax.add_patch(plt.Rectangle((i + 0.15 + spacing*i+0.2, values_ground5[i] - errors_ground5[i]), 0.15, 2*errors_ground5[i], color=color, alpha=0.7, hatch=hatches[i%len(hatches)]))
    ax.plot([i + spacing*i, i + spacing*i], [values_ground5[i] - errors_ground5[i], values_ground5[i] + errors_ground5[i]], color='blue', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M4_first
for i in range(len(values_first3)):
    if i < 6:
        color = avg_8_10[i]
    else:
        color = avg_8_10[i-6]
    ax.add_patch(plt.Rectangle((i + 0.15 + spacing*i+0.2, values_first5[i] - errors_first5[i]), 0.15, 2*errors_first5[i], color=color, alpha=0.7, hatch=hatches[i%len(hatches)]))
    ax.plot([i + spacing*i, i + spacing*i], [values_first5[i] - errors_first5[i], values_first5[i] + errors_first5[i]], color='green', alpha=0.0)

# Plot rectangles for each value and error_diff with a bit of spacing for M4_second
for i in range(len(values_second3)):
    if i < 6:
        color = avg_8_10[i]
    else:
        color = avg_8_10[i-6]
    ax.add_patch(plt.Rectangle((i + 0.15 + spacing*i+0.2, values_second5[i] - errors_second5[i]), 0.15, 2*errors_second5[i], color=color, alpha=0.7, hatch=hatches[i%len(hatches)]))
    ax.plot([i + spacing*i, i + spacing*i], [values_second5[i] - errors_second5[i], values_second5[i] + errors_second5[i]], color='black', alpha=0.0)

# Set labels and title
#ax.set_xlabel('Channels', y=-0.80, fontsize=14)
ax.set_ylabel('Energy levels $[w_0^{-1}]$', fontsize=14)
#ax.set_title('Mesonic spectrum, $am^{\\rm f}_0 = -0.71$, $am^{\\rm as}_0 = -1.01$', fontsize=16)
ax.set_title('Meson spectrum', fontsize=16)

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



# Plot transparent colorblocks to delimit xticks
for i in range(len(x_labels)):
    ax.add_patch(plt.Rectangle((i + spacing*i - 0.7 +0.2, ax.get_ylim()[0]), 1.0, ax.get_ylim()[1] - ax.get_ylim()[0], color='gray', alpha=0.08))
    
   
# Create custom legend handles with color patches for each group
legend_handles = [
    mpatches.Patch(color=avg_1_4[0], alpha=0.7),
    mpatches.Patch(color=avg_1_4[1], alpha=0.7),
    mpatches.Patch(color=avg_1_4[2], alpha=0.7),
    mpatches.Patch(color=avg_1_4[3], alpha=0.7),
    mpatches.Patch(color=avg_1_4[4], alpha=0.7),
    mpatches.Patch(color=avg_1_4[5], alpha=0.7),
    mpatches.Patch(color=avg_2_5[0], alpha=0.7),
    mpatches.Patch(color=avg_2_5[1], alpha=0.7),
    mpatches.Patch(color=avg_2_5[2], alpha=0.7),
    mpatches.Patch(color=avg_2_5[3], alpha=0.7),
    mpatches.Patch(color=avg_2_5[4], alpha=0.7),
    mpatches.Patch(color=avg_2_5[5], alpha=0.7),
    mpatches.Patch(color=avg_3_6[0], alpha=0.7),
    mpatches.Patch(color=avg_3_6[1], alpha=0.7),
    mpatches.Patch(color=avg_3_6[2], alpha=0.7),
    mpatches.Patch(color=avg_3_6[3], alpha=0.7),
    mpatches.Patch(color=avg_3_6[4], alpha=0.7),
    mpatches.Patch(color=avg_3_6[5], alpha=0.7),
    mpatches.Patch(color=avg_7_9[0], alpha=0.7),
    mpatches.Patch(color=avg_7_9[1], alpha=0.7),
    mpatches.Patch(color=avg_7_9[2], alpha=0.7),
    mpatches.Patch(color=avg_7_9[3], alpha=0.7),
    mpatches.Patch(color=avg_7_9[4], alpha=0.7),
    mpatches.Patch(color=avg_7_9[5], alpha=0.7),
    mpatches.Patch(color=avg_8_10[0], alpha=0.7),
    mpatches.Patch(color=avg_8_10[1], alpha=0.7),
    mpatches.Patch(color=avg_8_10[2], alpha=0.7),
    mpatches.Patch(color=avg_8_10[3], alpha=0.7),
    mpatches.Patch(color=avg_8_10[4], alpha=0.7),
    mpatches.Patch(color=avg_8_10[5], alpha=0.7),
]

# Create custom legend labels for each group
legend_labels = ['M4', 'M2', 'M3'] * 2  # Repeated for each color pair

# Plot the legend with custom handles and labels
#plt.legend(handles=legend_handles, labels=legend_labels, handlelength=1.5, handletextpad=3.4)
plt.legend([(legend_handles[0], legend_handles[1], legend_handles[2], legend_handles[3], legend_handles[4], legend_handles[5]), (legend_handles[6], legend_handles[7], legend_handles[8], legend_handles[9], legend_handles[10], legend_handles[11]), (legend_handles[12],legend_handles[13],legend_handles[14],legend_handles[15],legend_handles[16],legend_handles[17]), (legend_handles[18], legend_handles[19],legend_handles[20],legend_handles[21],legend_handles[22],legend_handles[23]), (legend_handles[24], legend_handles[25],legend_handles[26],legend_handles[27],legend_handles[28],legend_handles[29])], ['M1', 'M2', 'M3', 'M4', 'M5'],
               handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=4.5)



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
plt.ylim(0.79,3.70)

plt.tight_layout()
plt.savefig(output_file, format='pdf')
#plt.show()

