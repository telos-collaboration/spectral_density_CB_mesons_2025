import matplotlib.pyplot as plt
import numpy as np

plt.style.use("paperdraft.mplstyle")


plt.figure(figsize=(4.5, 8.0))

# Read data from the first file
filename1 = 'varying_sigma_s0p3.txt'
data1 = np.loadtxt(filename1)

# Read data from the second file
filename2 = 'varying_sigma_s0p1.txt'
data2 = np.loadtxt(filename2)

# Read data from the third file
filename3 = 'varying_sigma_s0p7.txt'
data3 = np.loadtxt(filename3)

mpi = 0.33

# Extract columns for the first file
x1 = data1[:, 0] / mpi
y1 = data1[:, 2]
y_error1 = data1[:, -1]

# Extract columns for the second file
x2 = data2[:, 0] / mpi
y2 = data2[:, 2]
y_error2 = data2[:, -1]

# Extract columns for the third file
x3 = data3[:, 0] / mpi
y3 = data3[:, 2]
y_error3 = data3[:, -1]

# Calculate scaling factors based on the integrals
integral1 = np.trapz(y1, x=x1)
integral2 = np.trapz(y2, x=x2)
integral3 = np.trapz(y3, x=x3)

scaling_factor2 = integral1 / integral2
scaling_factor3 = integral1 / integral3

# Scale the y-height and errors of the second dataset
y2_scaled = y2 * scaling_factor2 * 0.6
y_error2_scaled = y_error2 * scaling_factor2 * 0.6

# Scale the y-height and errors of the third dataset
y3_scaled = y3 * scaling_factor3
y_error3_scaled = y_error3 * scaling_factor3

# Create subplots with three panels
fig, axes = plt.subplots(3, 1, figsize=(6.5, 7))


# Plot the scaled data with error bars and lines for the second file
axes[0].errorbar(x2, y2_scaled, yerr=y_error2_scaled, fmt='o', label='$\sigma = 0.1m_0$', elinewidth=1.4, capsize=4, markersize=5)
axes[0].plot(x2, y2_scaled, linestyle='-', color='C0', linewidth=1.4)
axes[0].set_ylabel(' ', fontsize=14)
axes[0].legend()

# Plot the data with error bars and lines for the first file
axes[1].errorbar(x1, y1, yerr=y_error1*0.8, fmt='o', label='$\sigma = 0.3m_0$', color='C1', elinewidth=1.4, capsize=4, markersize=5)
axes[1].plot(x1, y1, linestyle='-', color='C1', linewidth=1.4)
axes[1].set_ylabel('$\hat{\\rho}_\sigma(\omega)$', fontsize=16)
axes[1].legend()

axes[1].yaxis.set_label_coords(-0.08, 0.5)  # Adjust the vertical position by modifying the second argument

# Plot the scaled data with error bars and lines for the third file
axes[2].errorbar(x3, y3_scaled, yerr=y_error3_scaled, fmt='o', label='$\sigma = 0.7m_0$',color='C2', elinewidth=1.4, capsize=4, markersize=5)
axes[2].plot(x3, y3_scaled, linestyle='-', color='C2', linewidth=1.4)
axes[2].set_xlabel('$\omega/m_0$', fontsize=16)
axes[2].set_ylabel(' ', fontsize=14)
axes[2].legend()

axes[2].xaxis.set_label_coords(0.5, -0.23)  # Adjust the vertical position by modifying the second argument


axes[0].set_xticklabels([])  # Remove ticks and labels from the first subplot
axes[1].set_xticklabels([])  # Remove ticks and labels from the second subplot
axes[0].grid(linestyle='--')
axes[1].grid(linestyle='--')
axes[2].grid(linestyle='--')

# Add title to the entire figure
fig.suptitle('$m_0 = 0.33,\,\, m_0^{*} = 2m_0$', fontsize=15, y = 0.93)

plt.savefig("../../../plots/sigma_variation.pdf", dpi=300, bbox_inches='tight')

# Show the plot
#plt.show()
