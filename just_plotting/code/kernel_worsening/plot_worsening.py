import numpy as np
import matplotlib.pyplot as plt
import math

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

plt.style.use("paperdraft.mplstyle")

# Assuming you have the necessary data in Kernel1.txt and Kernel2.txt
# Load data from files
data_kernel1 = np.loadtxt('kernel_0.4211936842105263.txt')
data_kernel2 = np.loadtxt('kernel_0.6274926315789473.txt')
data_kernel3 = np.loadtxt('kernel_0.782216842105263.txt')

# Extract energies and y-values from the data
energies = np.linspace(0.1 * 0.4083, 3.0 * 0.4083, len(data_kernel1))
x_to_plot_kernel1 = data_kernel1[:, 0]
y_to_plot_kernel1 = data_kernel1[:, 1]
x_to_plot_kernel2 = data_kernel2[:, 0]
y_to_plot_kernel2 = data_kernel2[:, 1]
x_to_plot_kernel3 = data_kernel3[:, 0]
y_to_plot_kernel3 = data_kernel3[:, 1]

# Other parameters
massNorm = 0.4083
omega =  0.4211936842105263
omega2 =  0.6274926315789473
omega3 =  0.782216842105263

def halfnorm_fp(e, s):  # int_0^inf dE exp{(-e-e0)^2/2s^2}
    res_ = math.erf(e / (np.sqrt(2) * s))
    res_ += 1
    res_ *= s
    res_ *= np.sqrt(np.pi / 2)
    return res_

def gauss_fp(x, x0, sigma, norm="Full"):
    if sigma == 0:
        return kronecker_fp(x, x0)
    if norm == "Full" or norm == "full":
        return (np.exp(-((x - x0) ** 2) / (2 * (sigma**2)))) / (
            sigma * math.sqrt(2 * math.pi)
        )
    if norm == "None" or norm == "none":
        return np.exp(-((x - x0) ** 2) / (2 * sigma**2))
    if norm == "Half" or norm == "half":
        return (np.exp(-((x - x0) ** 2) / (2 * sigma**2))) / halfnorm_fp(x0, sigma)

# Create subplots
fig, axs = plt.subplots(3, 1, figsize=(6.9, 5.5))

# Plot the y_to_plot curves
y_to_plot1 = gauss_fp(energies, omega, 0.28*massNorm, norm="half")
y_to_plot2 = gauss_fp(energies, omega2, 0.28*massNorm, norm="half")
y_to_plot3 = gauss_fp(energies, omega3, 0.28*massNorm, norm="half")


axs[0].plot(
    energies / massNorm,
    y_to_plot1,
    ls='--',
    label='$\Delta_\sigma$ at $\omega / m_{\mathrm{V}} \\approx 1.0$',
    color='red',
    linewidth=1.8,
)
axs[1].plot(
    energies / massNorm,
    y_to_plot2,
    ls='-.',
    label='$\Delta_\sigma$ at $\omega / m_{\mathrm{V}} \\approx 1.5$',
    color='red',
    linewidth=1.8,
)
axs[2].plot(
    energies / massNorm,
    y_to_plot3,
    ls=':',
    label='$\Delta_\sigma$ at $\omega / m_{\mathrm{V}} \\approx 1.9$',
    color='red',
    linewidth=1.8,
)



# Plot Kernel1
axs[0].plot(
    x_to_plot_kernel1,
    y_to_plot_kernel1,
    marker="o",
    markersize=6.8,
    ls="--",
    label='$\\bar{\Delta}_\sigma$',
    linewidth=1.5,
    color='black',
    markerfacecolor=CB_color_cycle[0],  # Replace with the actual color
)


# Plot Kernel2
axs[1].plot(
    x_to_plot_kernel2,
    y_to_plot_kernel2,
    marker="o",
    markersize=6.8,
    ls="--",
    label='$\\bar{\Delta}_\sigma$',
    linewidth=1.5,
    color='black',
    markerfacecolor=CB_color_cycle[1],  # Replace with the actual color
)


# Plot Kernel3
axs[2].plot(
    x_to_plot_kernel3,
    y_to_plot_kernel3,
    marker="o",
    markersize=6.8,
    ls="--",
    label='$\\bar{\Delta}_\sigma$',
    linewidth=1.5,
    color='black',
    markerfacecolor=CB_color_cycle[2],  # Replace with the actual color
)



# Add legend and labels
axs[2].set_xlabel('$E/m_{\mathrm{V}}$', fontsize=14)
plt.tight_layout()
axs[0].legend()
axs[1].legend()
axs[2].legend()

axs[0].grid(linestyle='--')
axs[1].grid(linestyle='--')
axs[2].grid(linestyle='--')

plt.savefig(
    '../../../plots/kernel_worsening.pdf',
    dpi=130, bbox_inches='tight'
)

# Show the plot
#plt.show()

