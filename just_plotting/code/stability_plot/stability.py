import os
import matplotlib.pyplot as plt
import numpy as np

def plot_stability_multiple_alpha(log_files, save_plot=True, n_alphas=2):
    plt.style.use("paperdraft.mplstyle")
    fig, ax = plt.subplots(2, 1, figsize=(6, 8))

    plot_markers = ['o', 's', 'D']
    CB_colors = plt.get_cmap('tab10').colors[:10]

    alpha_values = [0.00, 1.00, 1.99]  # Manually set alpha values

    for i, log_file in enumerate(log_files[:n_alphas]):
        data = np.loadtxt(log_file, usecols=(1, 2, 4, 5))

        ax[0].errorbar(
            x=data[:, 0],
            y=data[:, 1],
            yerr=data[:, 2],
            marker=plot_markers[i],
            markersize=3.8,
            elinewidth=1.4,
            capthick=1.4,
            capsize=3,
            ls="",
            label=r"$\alpha = {:1.2f}$".format(alpha_values[i]),
            color=CB_colors[i],
        )

        if i == 0:
            ax[1].errorbar(
                x=data[:, 3],
                y=data[:, 1],
                yerr=data[:, 2],
                marker=plot_markers[i],
                markersize=3.8,
                capthick=1.4,
                elinewidth=1.4,
                capsize=3,
                ls="",
                label=r"$\alpha = {:1.2f}$".format(alpha_values[i]),
                color=CB_colors[i],
            )

            # Add color block for the sixth smallest lambda value
            min_lambda_index = np.argsort(data[:, 3])[5]
            mean_lambda = data[min_lambda_index, 3]
            width_lambda = data[min_lambda_index, 2]

            ax[0].axhspan(
                ymin=data[min_lambda_index, 1] - data[min_lambda_index, 2],
                ymax=data[min_lambda_index, 1] + data[min_lambda_index, 2],
                alpha=0.3,
                color=CB_colors[4],
            )
            
            ax[1].axhspan(
                ymin=data[min_lambda_index, 1] - data[min_lambda_index, 2],
                ymax=data[min_lambda_index, 1] + data[min_lambda_index, 2],
                alpha=0.3,
                color=CB_colors[4],
            )

    ax[0].set_xlabel(r"$\lambda$", fontsize=14)
    ax[0].set_ylabel(r"$\hat{\rho}_\sigma (\omega)$", fontsize=14)
    ax[0].legend(prop={"size": 14, "family": "Helvetica"})
    ax[0].set_xscale('log')
    ax[0].grid()

    ax[1].set_xlabel(r"$A[\vec{g}] / A_0$", fontsize=14)
    ax[1].set_ylabel(r"$\hat{\rho}_\sigma (\omega)$", fontsize=14)
    ax[1].legend(prop={"size": 14, "family": "Helvetica"})
    ax[1].set_xscale('log')
    ax[1].grid()

    plt.tight_layout()
    plt.suptitle('$\omega/m_{\mathrm{V}} = 1.25$ $\sigma = 0.30m_{\mathrm{V}}$',fontsize=14, y=0.99)
    if save_plot:
        plt.savefig("../../../plots/LambdaScan.pdf", dpi=300)
    
    #plt.show()

if __name__ == "__main__":
    log_files = ["InverseProblemLOG_AlphaA.log", "InverseProblemLOG_AlphaB.log", "InverseProblemLOG_AlphaC.log"]
    plot_stability_multiple_alpha(log_files, save_plot=True, n_alphas=3)
