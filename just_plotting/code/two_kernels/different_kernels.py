import numpy as np
import matplotlib.pyplot as plt
import math

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']


factor = 0.6


plt.style.use("paperdraft.mplstyle")
plt.figure(figsize=(5, 3.5))
# Assuming you have the necessary data in Kernel1.txt and Kernel2.txt
# Load data from files
data_kernel1 = np.loadtxt('gauss_kernel_0.2285857894736842.txt')
data_kernel2 = np.loadtxt('cauchy_kernel_0.2285857894736842.txt')

# Extract energies and y-values from the data
energies = np.linspace(0.05 * 0.4054, 2.65 * 0.4054 - 0.07, len(data_kernel1)*10)
x_to_plot_kernel1 =data_kernel1[:, 0]
y_to_plot_kernel1 = data_kernel1[:, 1]
x_to_plot_kernel2 =data_kernel2[:, 0]
y_to_plot_kernel2 = data_kernel2[:, 1]

# Other parameters

massNorm = 0.4054
omega =  0.2285857894736842

def cauchy(x, mu1, sigma):
    return sigma / (sigma**2 + (mu1 - x)**2)


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
    x_to_plot_kernel1,
    y_to_plot_kernel1,
    marker="o",
    markersize=6.8,
    ls="--",
    linewidth=1.6,
    label="Gauss Ker. at $\omega / m_\mathrm{V} = 0.6$",
    color='black',
    markerfacecolor=CB_color_cycle[1],  # Replace with the actual color
)

# Plot Kernel2
plt.plot(
    x_to_plot_kernel2,
    y_to_plot_kernel2*factor,
    marker="o",
    markersize=6.8,
    ls="--",
    linewidth=1.6,
    label="Cauchy Ker. at $\omega / m_\mathrm{V} = 0.6$",
    color='black',
    markerfacecolor=CB_color_cycle[0],  # Replace with the actual color
)

plt.title("$\quad \\alpha = 0.0 \quad \sigma = 0.40m_{\mathrm{V}}$", fontsize=13)
y_to_plot = gauss_fp(energies, omega, 0.40*massNorm, norm="half")
y_to_plot2 = cauchy(energies, omega, 0.40*massNorm)

plt.plot(
            energies / massNorm,
            y_to_plot,
            ls = '-',
            label = 'Gaussian target',
            color='blue',
            linewidth=1.8,
        )

plt.plot(
            energies / massNorm,
            y_to_plot2*factor,
            ls = '-',
            label = 'Cauchy target',
            color='red',
            linewidth=1.8,
        )


# Add legend and labels
plt.legend(fontsize=12,frameon=False)
plt.grid(linestyle='--')
plt.xlabel('$E/m_{\mathrm{V}}$', fontsize=14)
plt.tight_layout()

plt.savefig(
            '../../../plots/two_kernels_comparisons.pdf',
            dpi=130, bbox_inches='tight'
        )

# Show the plot
#plt.show()

