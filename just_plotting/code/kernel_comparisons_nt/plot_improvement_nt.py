import numpy as np
import matplotlib.pyplot as plt
import math

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']


plt.style.use("paperdraft.mplstyle")
plt.figure(figsize=(6, 4.0))
# Assuming you have the necessary data in Kernel1.txt and Kernel2.txt
# Load data from files
data_kernel1 = np.loadtxt('Kernel1.txt')
data_kernel2 = np.loadtxt('Kernel2.txt')
data_kernel3 = np.loadtxt('Kernel3.txt')

# Extract energies and y-values from the data
energies = np.linspace(0.05 * 0.4054, 4.0 * 0.4054, len(data_kernel1))
y_to_plot_kernel1 = data_kernel1[:, 1]
y_to_plot_kernel2 = data_kernel2[:, 1]
y_to_plot_kernel3 = data_kernel3[:, 1]

# Other parameters

massNorm = 0.4054
omega =  1.8*massNorm


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




# Plot Kernel1
plt.plot(
    energies / massNorm,
    y_to_plot_kernel1,
    marker="o",
    markersize=6.8,
    ls="--",
    linewidth=1.6,
    label="Gauss Ker., $N_t = 48$".format(omega / massNorm),
    color='black',
    markerfacecolor=CB_color_cycle[0],  # Replace with the actual color
)

# Plot Kernel2
plt.plot(
    energies / massNorm,
    y_to_plot_kernel2,
    marker="o",
    markersize=6.8,
    ls="--",
    linewidth=1.6,
    label="Gauss Ker., $N_t = 64$".format(omega / massNorm),
    color='black',
    markerfacecolor=CB_color_cycle[1],  # Replace with the actual color
)


# Plot Kernel3
plt.plot(
    energies / massNorm,
    y_to_plot_kernel3,
    marker="o",
    markersize=6.8,
    ls="--",
    linewidth=1.6,
    label="Gauss Ker., $N_t = 96$".format(omega / massNorm),
    color='black',
    markerfacecolor=CB_color_cycle[2],  # Replace with the actual color
)


plt.title("$\omega/m_V$ = " + "{:2.1f}".format(omega / massNorm) + "$\quad \\alpha = 0.0 \quad \sigma = 0.36m_V $", fontsize=14)
y_to_plot = gauss_fp(energies, omega, 0.36*massNorm, norm="half")

plt.plot(
            energies / massNorm,
            y_to_plot,
            ls = '-',
            label = 'Target',
            color='red',
            linewidth=1.8,
        )


# Add legend and labels
plt.legend(fontsize=11,frameon=False)
plt.xlabel('$E/m_V$', fontsize=14)
plt.tight_layout()

plt.savefig(
            '../../../plots/kernel_comparison.pdf',
            dpi=130, bbox_inches='tight'
        )

# Show the plot
#plt.show()

