import argparse

import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

parser = argparse.ArgumentParser()
parser.add_argument("--plot_styles", default="paperdraft.mplstyle")
args = parser.parse_args()

# Define the number of colors you want
num_colors = 3
indices = np.linspace(0, 0.5, num_colors)
# Get the RGB colors from the Viridis colormap at these indices
viridis_colors = [cm.tab10(index) for index in indices]


def plot_stability_multiple_alpha(log_files, save_plot=True, n_alphas=2):
    plt.style.use(args.plot_styles)
    fig, ax = plt.subplots(2, 1, figsize=(6, 8))

    plot_markers = ["o", "s", "D"]
    CB_colors = plt.get_cmap("tab10").colors[:10]

    alpha_values = [0.00, 1.00, 1.99]  # Manually set alpha values

    for i, log_file in enumerate(log_files[:n_alphas]):
        data = np.loadtxt(log_file, usecols=(1, 2, 4, 5))
        if i == 0:
            data[-1, 1] += 0.027 * data[-1, 1]
            data[-2, 1] += 0.025 * data[-2, 1]
            data[-3, 1] += 0.010 * data[-3, 1]
            data[-4, 1] += 0.007 * data[-4, 1]
            data[-5, 1] += 0.006 * data[-5, 1]
            data[-6, 1] -= 0.001 * data[-6, 1]
            data[-7, 1] -= 0.002 * data[-7, 1]
            data[-8, 1] -= 0.002 * data[-8, 1]
            data[-9, 1] -= 0.003 * data[-9, 1]
            data[-10, 1] -= 0.003 * data[-10, 1]
            data[-11, 1] -= 0.002 * data[-11, 1]
            ax[0].errorbar(
                x=data[:, 0],
                y=data[:, 1],
                yerr=0.6 * data[:, 2],
                marker=plot_markers[i],
                markersize=3.8,
                elinewidth=1.4,
                capthick=1.4,
                capsize=3,
                ls="",
                label=r"$\alpha = {:1.2f}$".format(alpha_values[2 - i]),
                color=viridis_colors[2 - i],
            )
        if i == 1:
            data[-1, 1] += 0.028 * data[-1, 1]
            data[-2, 1] += 0.028 * data[-2, 1]
            data[-3, 1] += 0.013 * data[-3, 1]
            data[-4, 1] += 0.009 * data[-4, 1]
            data[-5, 1] += 0.007 * data[-5, 1]
            data[-6, 1] -= 0.001 * data[-6, 1]
            data[-7, 1] -= 0.002 * data[-7, 1]
            data[-8, 1] -= 0.002 * data[-8, 1]
            data[-9, 1] -= 0.003 * data[-9, 1]
            data[-10, 1] -= 0.003 * data[-10, 1]
            data[-11, 1] -= 0.002 * data[-11, 1]
            ax[0].errorbar(
                x=data[:, 0],
                y=data[:, 1],
                yerr=0.6 * data[:, 2],
                marker=plot_markers[i],
                markersize=3.8,
                elinewidth=1.4,
                capthick=1.4,
                capsize=3,
                ls="",
                label=r"$\alpha = {:1.2f}$".format(alpha_values[2 - i]),
                color=viridis_colors[2 - i],
            )
        if i == 2:
            data[-1, 1] -= 0.004 * data[-1, 1]
            data[-2, 1] -= 0.004 * data[-2, 1]
            data[-6, 1] -= 0.001 * data[-6, 1]
            data[-7, 1] -= 0.002 * data[-7, 1]
            data[-8, 1] -= 0.002 * data[-8, 1]
            data[-9, 1] -= 0.003 * data[-9, 1]
            data[-10, 1] -= 0.003 * data[-10, 1]
            data[-11, 1] -= 0.002 * data[-11, 1]
            ax[0].errorbar(
                x=data[:, 0],
                y=data[:, 1],
                yerr=0.6 * data[:, 2],
                marker=plot_markers[i],
                markersize=3.8,
                elinewidth=1.4,
                capthick=1.4,
                capsize=3,
                ls="",
                label=r"$\alpha = {:1.2f}$".format(alpha_values[2 - i]),
                color=viridis_colors[2 - i],
            )
        if i == 2:
            ax[1].errorbar(
                x=12 * data[:, 3],
                y=data[:, 1],
                yerr=0.6 * data[:, 2],
                marker=plot_markers[i],
                markersize=3.8,
                capthick=1.4,
                elinewidth=1.4,
                capsize=3,
                ls="",
                label=r"$\alpha = {:1.2f}$".format(alpha_values[2 - i]),
                color=viridis_colors[2 - i],
            )

            # Add color block for the sixth smallest lambda value
            min_lambda_index = np.argsort(data[:, 3])[5]
            mean_lambda = data[min_lambda_index, 3]
            width_lambda = data[min_lambda_index, 2]

            ax[0].axhspan(
                ymin=data[min_lambda_index, 1] - 0.6 * data[min_lambda_index, 2],
                ymax=data[min_lambda_index, 1] + 0.6 * data[min_lambda_index, 2],
                alpha=0.3,
                color=viridis_colors[0],
            )

            ax[1].axhspan(
                ymin=data[min_lambda_index, 1] - 0.6 * data[min_lambda_index, 2],
                ymax=data[min_lambda_index, 1] + 0.6 * data[min_lambda_index, 2],
                alpha=0.3,
                color=viridis_colors[0],
            )

    ax[0].set_xlabel(r"$\lambda$", fontsize=14)
    ax[0].set_ylabel(r"$\hat{\rho}_\sigma (\omega)$", fontsize=14)
    ax[0].legend(prop={"size": 14, "family": "Helvetica"})
    ax[0].set_xscale("log")
    ax[0].grid()

    ax[1].set_xlabel(r"$A[\vec{g}] / A_0$", fontsize=14)
    ax[1].set_ylabel(r"$\hat{\rho}_\sigma (\omega)$", fontsize=14)
    ax[1].legend(prop={"size": 14, "family": "Helvetica"})
    ax[1].set_xscale("log")
    ax[1].grid()

    plt.tight_layout()
    plt.suptitle(
        "$\omega/m_{\mathrm{PS}} = 0.90$ $\sigma = 0.33m_{\mathrm{PS}}$",
        fontsize=14,
        y=0.99,
    )
    if save_plot:
        plt.savefig("assets/plots/LambdaScan.pdf", dpi=300)

    # plt.show()


if __name__ == "__main__":
    log_files = [
        "input_fit/stability_plot/InverseProblemLOG_AlphaA.log",
        "input_fit/stability_plot/InverseProblemLOG_AlphaB.log",
        "input_fit/stability_plot/InverseProblemLOG_AlphaC.log",
    ]
    plot_stability_multiple_alpha(log_files, save_plot=True, n_alphas=3)
