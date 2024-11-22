import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.legend_handler import HandlerTuple
import matplotlib.patches as mpatches

# Define constants and configuration
w0_values = {'M1': 2.5210, 'M2': 2.5290, 'M3': 2.5237, 'M4': 2.3664, 'M5': 2.6927}
files = [
    ('M1', ['M1_ground.txt', 'M1_first.txt']),
    ('M2', ['M2_ground.txt', 'M2_first.txt', 'M2_second.txt']),
    ('M3', ['M3_ground.txt', 'M3_first.txt', 'M3_second.txt']),
    ('M4', ['M4_ground.txt', 'M4_first.txt', 'M4_second.txt']),
    ('M5', ['M5_ground.txt', 'M5_first.txt', 'M5_second.txt']),
    # CB channels
    ('M1', ['CB_M1_ground.txt', 'CB_M1_first.txt', 'CB_M1_second.txt']),
    ('M2', ['CB_M2_ground.txt', 'CB_M2_first.txt', 'CB_M2_second.txt']),
    ('M3', ['CB_M3_ground.txt', 'CB_M3_first.txt', 'CB_M3_second.txt']),
    ('M4', ['CB_M4_ground.txt', 'CB_M4_first.txt', 'CB_M4_second.txt']),
    ('M5', ['CB_M5_ground.txt', 'CB_M5_first.txt', 'CB_M5_second.txt'])
]

plt.style.use("paperdraft.mplstyle")
output_file = '../../../plots/final_spectrum_MN3.pdf'
spacing = 0.
colors = [cm.tab10(i) for i in np.linspace(0, 0.5, 6)]

def process_file(file_name, w0):
    values, errors = [], []
    with open(file_name, 'r') as file:
        for line in file:
            data = list(map(float, line.strip().split()))
            measurements, errors_data = data[::2], data[1::2]
            min_err_idx = errors_data.index(min(errors_data))
            stat = errors_data[min_err_idx]
            sys = max(abs(measurements[i] - measurements[j]) for i in range(len(measurements)) for j in range(i+1, len(measurements)))
            value = measurements[min_err_idx] * w0
            total_err = math.sqrt(stat**2 + sys**2) * w0
            values.append(value)
            errors.append(total_err)
    return values, errors

# Plotting setup
fig, ax = plt.subplots(figsize=(16, 9))

# Define x_labels: Place 'η'' first, followed by the rest of the labels
x_labels = ['PS', 'V', 'T', 'AV', 'AT', 'S', 'ps', 'v', 't', 'av', 'at', 's', 
            '$\Lambda^{+}_{\\rm CB}$', '$\Sigma^{+}_{\\rm CB}$', '$\Sigma^{*, +}_{\\rm CB}$', 
            '$\Lambda^{-}_{\\rm CB}$', '$\Sigma^{-}_{\\rm CB}$', '$\Sigma^{*, -}_{\\rm CB}$']

hatches = [r'//////', '++++']  # Define different hatch patterns for M4 and M5

# Offsets and alpha values for ensembles
offsets = {'M1': -0.05, 'M2': 0.0, 'M3': 0.05, 'M4': -0.25, 'M5': 0.50}
alpha_values = {'M1': 0.4, 'M2': 0.6, 'M3': 0.85, 'M4': 0.7, 'M5': 0.7}

# Define the new color scheme for CB channels
cb_colors = [cm.RdBu(i) for i in np.linspace(0, 0.3, 6)]  # Use a range from the RdBu colormap

# Define hatch colors based on the same colormap as the bars
hatch_colors = [cm.tab10(i) for i in np.linspace(0, 0.5, 6)]


# Plot the rest of the channels
for ensemble, file_list in files:
    w0 = w0_values[ensemble]
    ensemble_offset = offsets[ensemble]

    for idx, file_name in enumerate(file_list):
        values, errors = process_file(file_name, w0)
        for i, (val, err) in enumerate(zip(values, errors)):
            if 'CB' in file_name:
                x_pos = 13.1 + i + ensemble_offset
                color = cb_colors[i % 6]
            else:
                x_pos = 0.3 + i + ensemble_offset
                color = colors[i % 6]
                if i >= 6:
                    x_pos += 0.4

            hatch = None
            hatch_color = None
            alpha = alpha_values.get(ensemble, 0.6)

            if ensemble == 'M4':
                hatch = hatches[0]
                hatch_color = colors[i % 6]
                color = 'none'
                if 'CB' in file_name:
                    hatch_color = cb_colors[i % 6]
            elif ensemble == 'M5':
                hatch = hatches[1]
                hatch_color = colors[i % 6]
                color = 'none'
                if 'CB' in file_name:
                    hatch_color = cb_colors[i % 6]

            rect = plt.Rectangle((x_pos - 0.2, val - err), 0.15 if ensemble in ['M4', 'M5'] else 0.4, 2 * err, color=color, alpha=alpha, hatch=hatch)
            if hatch_color:
                rect.set_edgecolor(hatch_color)
            ax.add_patch(rect)

# Define sector boundaries and labels
sector_boundaries = [6.0, 12.4, 18.8]
sector_labels = ['Fundamental', 'Antisymmetric', 'Chimera baryons']

# Add vertical lines for sector separation
for boundary in sector_boundaries:
    ax.axvline(boundary, color='black', linestyle='-', alpha=0.7)

# Optionally shade sectors to distinguish regions
ax.axvspan(-0.5, 7.5, color='none')  # First sector
ax.axvspan(0.5, 7.5, color='none')  # Second sector
ax.axvspan(7.5, 13.5, color='none')  # Third sector

# Add labels above the sectors
for i, label in enumerate(sector_labels):
    x_pos = (sector_boundaries[i - 1] + sector_boundaries[i]) / 2 if i > 0 else (sector_boundaries[i] - 0.5) / 2
    ax.text(x_pos, 4.1, label, fontsize=16, ha='center', va='center', transform=ax.transData)

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

ciao = [0.3, 1.3, 2.3, 3.3, 4.3, 5.3, 6.6, 7.6, 8.6, 9.6, 10.6, 11.6, 13.1, 14.1, 15.1, 16.1, 17.1, 18.1]
# Adjust x-ticks and labels
ax.set_xticks(ciao)
ax.set_xticklabels(x_labels, ha='right', fontsize=14, rotation=45)

# Set y-label and axis limits
ax.set_ylabel('$\hat{m}$', fontsize=16)
ax.set_ylim(0.5, 4.3)
ax.set_xlim(-0.5, 18.8)
# Save and show the plot
plt.tight_layout()
plt.savefig(output_file)
#plt.show()

