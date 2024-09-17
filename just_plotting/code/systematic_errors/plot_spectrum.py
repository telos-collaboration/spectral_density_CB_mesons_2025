import numpy as np
import matplotlib.pyplot as plt

# Activating text rendering by LaTeX
plt.style.use("paperdraft.mplstyle")

# Cool colors for the four points
cool_colors = ['#1f78b4', '#33a02c', '#6a3d9a', '#e31a1c']

# Marker styles for Cauchy points
cauchy_marker_styles = ['o', 's', 'D', '^', 'v', '<']

# Read data from the file
file_path = 'spectrum_ground_fund.txt'
data = np.loadtxt(file_path)

x = np.arange(1, 7)
x_labels = ['PS', 'V', 'T', 'AV', 'AT', 'S']

# Create subplots with 1 row and 2 columns
fig, axs = plt.subplots(1, 2, figsize=(9, 4))

# Define the horizontal offsets
offset_gauss_k = -0.20
offset_gauss_k1 = 0.20
offset_cauchy_k = -0.065
offset_cauchy_k1 = 0.065

legend_handles = []  # Store legend handles for the common legend

for i in range(3):
    # Extracting y values and errors for each row
    y1 = data[i, 0:2]
    y2 = data[i, 2:4]
    y3 = data[i, 4:6]
    y4 = data[i, 6:8]

    # Plot y1 and y2 in the first subplot with offset for Gauss peaks
    handle1 = axs[0].errorbar(x[i] + offset_gauss_k, y1[0], yerr=y1[1], fmt='o', markersize=5, capsize=3.5, linewidth=2, color=cool_colors[0], markeredgecolor='black', label='$k$ peaks - Gauss' if i == 0 else None)
    handle2 = axs[0].errorbar(x[i] + offset_gauss_k1, y2[0], yerr=y2[1], fmt='o', markersize=5, capsize=3.5, linewidth=2, color=cool_colors[1], markeredgecolor='black', label='$k+1$ peaks - Gauss' if i == 0 else None)
    
    # Plot y3 and y4 in the same subplot with offset for Cauchy peaks
    handle3 = axs[0].errorbar(x[i] + offset_cauchy_k, y3[0], yerr=y3[1], fmt=cauchy_marker_styles[2], markersize=5, capsize=3.5, linewidth=2, color=cool_colors[2], markeredgecolor='black', label='$k$ peaks - Cauchy' if i == 0 else None)
    handle4 = axs[0].errorbar(x[i] + offset_cauchy_k1, y4[0], yerr=y4[1], fmt=cauchy_marker_styles[2], markersize=5, capsize=3.5, linewidth=2, color=cool_colors[3], markeredgecolor='black', label='$k+1$ peaks - Cauchy' if i == 0 else None)

    # Store legend handles for the common legend
    if i == 0:
        legend_handles.extend([handle1, handle2, handle3, handle4])

    # Add shading between data points
    axs[0].fill_betweenx([0, 1], x[i] + 1.2*offset_gauss_k, x[i] + 1.2*offset_gauss_k1, color='gray', alpha=0.2)

axs[0].set_xticks(x[:3])
axs[0].set_xticklabels(x_labels[:3])
axs[0].grid(True, linestyle='--', alpha=0.7)

for i in range(3, 6):
    # Extracting y values and errors for each row
    y1 = data[i, 0:2]
    y2 = data[i, 2:4]
    y3 = data[i, 4:6]
    y4 = data[i, 6:8]

    # Plot y1 and y2 in the second subplot with offset for Gauss peaks
    axs[1].errorbar(x[i] + offset_gauss_k, y1[0], yerr=y1[1], fmt='o', markersize=5, capsize=3.5, linewidth=2, color=cool_colors[0], markeredgecolor='black', label='$k$ peaks - Gauss' if i == 3 else None)
    axs[1].errorbar(x[i] + offset_gauss_k1, y2[0], yerr=y2[1], fmt='o', markersize=5, capsize=3.5, linewidth=2, color=cool_colors[1], markeredgecolor='black', label='$k+1$ peaks - Gauss' if i == 3 else None)
    
    # Plot y3 and y4 in the same subplot with offset for Cauchy peaks
    axs[1].errorbar(x[i] + offset_cauchy_k, y3[0], yerr=y3[1], fmt=cauchy_marker_styles[2], markersize=5, capsize=3.5, linewidth=2, color=cool_colors[2], markeredgecolor='black', label='$k$ peaks - Cauchy' if i == 3 else None)
    axs[1].errorbar(x[i] + offset_cauchy_k1, y4[0], yerr=y4[1], fmt=cauchy_marker_styles[2], markersize=5, capsize=3.5, linewidth=2, color=cool_colors[3], markeredgecolor='black', label='$k+1$ peaks - Cauchy' if i == 3 else None)

    # Add shading between data points
    axs[1].fill_betweenx([0, 1], x[i] + 1.2*offset_gauss_k, x[i] + 1.2*offset_gauss_k1, color='gray', alpha=0.2)

axs[1].set_xticks(x[3:])
axs[1].set_xticklabels(x_labels[3:])
axs[0].set_ylabel('$aE_0$', fontsize=14, labelpad=12)
axs[1].grid(True, linestyle='--', alpha=0.7)

axs[0].set_ylim((0.357,0.413))
axs[1].set_ylim((0.503,0.564))

# Add a common legend below both subplots
fig.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, 0.03), frameon=False, fontsize=14, ncol=4)

# Add a main title for the entire figure
#plt.suptitle('$aE_0$ fits, fundamental sector', fontsize=14, y=0.915)

# Adjust layout
plt.tight_layout(rect=[0, 0, 1, 0.96])

# Save the figure in PDF format with dpi=300 and specified size
plt.savefig('../../../plots/nt64_spectrum_ground_fund.pdf', dpi=300, bbox_inches='tight')

# Show the plot
# plt.show()


