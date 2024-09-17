import matplotlib.pyplot as plt
import numpy as np

plt.style.use("paperdraft.mplstyle")

def plot_color_bands(ax, centers, widths, colors, labels, heights):
    for center, width, color, label, height in zip(centers, widths, colors, labels, heights):
        ax.bar(center, height, width=2*width, color=color, edgecolor='black', linewidth=1.5, alpha=0.65, label=label)

def main():
    # Read data from the TXT file
    file_path = 'improving_spectrum_nt.txt'
    data = np.loadtxt(file_path)

    num_rows, num_columns = data.shape
    
    # Create a single row of subplots
    fig, axs = plt.subplots(1, num_rows, figsize=(5 * num_rows, 3), sharey=True)  # Change the width to 6 inches

    for i in range(num_rows):
        # Extract data for each row
        row_data = data[i, :]

        # Divide the row into center, width, and color
        centers = row_data[0::2]
        widths = row_data[1::2]
        colors = ['blue', 'green', 'orange']  # Add more colors as needed
        labels = ['$N_t = 48$', '    $N_t = 64   $', '    $N_t = 96$']  # Add labels for each color band with added spaces
        heights = [1.0, 0.8, 0.5]  # Add heights for each color band (adjust as needed)

        # Plot color bands for each subplot
        plot_color_bands(axs[i], centers, widths, colors, labels, heights)
        axs[i].set_yticks([])  # Hide x-axis ticks
        axs[i].set_yticklabels([])  # Hide x-axis tick labels

    axs[0].set_xlabel('$aE_0$', fontsize=14)
    axs[1].set_xlabel('$aE_1$', fontsize=14)
    axs[0].set_ylim(0, 1)  # Set the y-axis limits for all subplots
    axs[0].set_xlim(0.360, 0.3735)
    axs[1].set_xlim(0.610, 0.700)
    
    # Add a common legend below the subplots at a specific height (e.g., -0.3) with increased spacing
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=len(colors), handlelength=1, handletextpad=0.2)

    # Adjust layout and display the plot
    plt.tight_layout()
    # plt.show()

    plt.savefig(
        '../../../plots/improving_spectrum_nt.pdf',
        dpi=130, bbox_inches='tight'
    )

if __name__ == "__main__":
    main()

