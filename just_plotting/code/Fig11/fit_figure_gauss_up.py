import numpy as np
import matplotlib.pyplot as plt

plt.style.use("paperdraft.mplstyle")
plt.figure(figsize=(7, 4.5))
# Load data from gauss_curves.txt
gauss_data = np.loadtxt("gauss_curves_up.txt")

# Extracting mean, amplitude, and error values for each gaussian curve
mean1, amplitude1, err_amplitude1 = gauss_data[0]
mean2, amplitude2, err_amplitude2 = gauss_data[1]

mpi = 0.4098

sigma = 0.30*mpi

# Define gaussian functions

def gaussian(x, mean, amplitude):
    return amplitude * np.exp(-((x - mean)**2) / (2*sigma**2))

# Generate x values for plotting
x_values = np.linspace(0, 2.15,1000)

# Plot gaussian curves
plt.plot(x_values / mpi, gaussian(x_values, mean1, amplitude1), color="chocolate", linewidth=1.6)
plt.fill_between(x_values / mpi, gaussian(x_values, mean1, amplitude1 - err_amplitude1),
                 gaussian(x_values, mean1, amplitude1 + err_amplitude1), alpha=0.2, color='chocolate')

plt.plot(x_values / mpi, gaussian(x_values, mean2, amplitude2), color="olive", linewidth=1.6)
plt.fill_between(x_values / mpi, gaussian(x_values, mean2, amplitude2 - err_amplitude2),
                 gaussian(x_values, mean2, amplitude2 + err_amplitude2), alpha=0.3, color='olive')
                 
               
plt.plot(x_values / mpi, gaussian(x_values, mean1, amplitude1) + gaussian(x_values, mean2, amplitude2), color="orange", linewidth=1.6)

upper_band = [0] * len(x_values)
lower_band = [0] * len(x_values)

for j in range(len(x_values)):
    lower_band[j] = gaussian(x_values[j], mean1, amplitude1 - err_amplitude1) + gaussian(x_values[j], mean2, amplitude2 - err_amplitude2)
    upper_band[j] = gaussian(x_values[j], mean1, amplitude1 + err_amplitude1) + gaussian(x_values[j], mean2, amplitude2 + err_amplitude2)


plt.fill_between(x_values / mpi, lower_band, upper_band, alpha=0.25, color='orange')

# Load data from sp_dens_datapoints.txt
data_points = np.loadtxt("sp_dens_datapoints_gauss_up.txt")

# Extracting energy and spdens values
energy, spdens, err_spdens = data_points[:,0], data_points[:,1], data_points[:,2]

# Plot data points
plt.errorbar(energy/ mpi, spdens, yerr=err_spdens, fmt='o', color='black', markersize=3.0, elinewidth=1.2,)

# Set labels and legend
plt.xlabel('$\omega / m_{\mathrm{V}}$', fontsize=14)
plt.ylabel('$\hat{\\rho}_{\sigma}(\omega)$', fontsize=14)
#plt.title('Gauss kernel', fontsize=13)
plt.xticks(np.arange(0.5, 2.01, 0.5))
plt.xlim(0.10, 2.50)
plt.grid(linestyle='--')

# Show plot
#plt.show()

plt.savefig(
            '../../../plots/fit_PAUL_gi_fund_N20_N20_double_2.pdf',
            dpi=130, bbox_inches='tight'
        )
