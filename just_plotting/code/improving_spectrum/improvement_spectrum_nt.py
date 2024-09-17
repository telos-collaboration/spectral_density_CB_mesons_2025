import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter

plt.style.use("paperdraft.mplstyle")

def plot_error_bars(ax, centers, errors, colors, labels, heights, errorbar_linewidth):
    for center, error, color, label, height in zip(centers, errors, colors, labels, heights):
        ax.errorbar(center, height, xerr=error, fmt='o', color=color, capsize=5, label=label, elinewidth=errorbar_linewidth)

def main():
    # Read data from the TXT file
    file_path = 'improving_spectrum_nt.txt'
    data = np.loadtxt(file_path)

    num_rows, num_columns = data.shape
    
    # Create a single row of subplots
    fig, axs = plt.subplots(1, num_rows, figsize=(4.6 * num_rows, 3), sharey=True)  # Change the width to 6 inches

    errorbar_linewidth = 1.05  # Define the width of error bars

    for i in range(num_rows):
        # Extract data for each row
        row_data = data[i, :]

        # Divide the row into center, error, and color
        centers = row_data[0::2]
        errors = row_data[1::2]
        colors = ['blue', 'blue', 'blue']  # Add more colors as needed
        labels = ['$N_t = 48$', '$N_t = 64$', '$N_t = 96$']  # Add labels for each color band
        heights = np.arange(1, num_columns + 1) * 0.1  # Add heights for each point with an offset

        # Plot points with error bars for each subplot
        plot_error_bars(axs[i], centers, errors, colors, labels, heights, errorbar_linewidth)
        axs[i].set_yticks([])  # Hide y-axis ticks
        axs[i].set_yticklabels([])  # Hide y-axis tick labels

    # Set y-axis labels for the first subplot
    tick_locations = heights[:3]   # Compute tick locations based on heights
    axs[0].set_yticks(tick_locations)
    axs[0].set_yticklabels(labels)

    # Specify formatting for x-axis tick labels
    for ax in axs:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    axs[0].set_xlabel('$aE_0$', fontsize=14)
    axs[1].set_xlabel('$aE_1$', fontsize=14)
    axs[0].set_ylim(0, num_columns * 0.07)  # Set the y-axis limits for all subplots with an additional buffer
    axs[0].set_xlim(0.360, 0.3735)
    axs[1].set_xlim(0.610, 0.700)

    # Adjust layout and display the plot
    plt.tight_layout()
    # plt.show()

    plt.savefig(
        '../../../plots/improving_spectrum_nt.pdf',
        dpi=130, bbox_inches='tight'
    )

if __name__ == "__main__":
    main()

