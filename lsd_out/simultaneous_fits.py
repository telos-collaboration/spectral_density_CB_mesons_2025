import os
import re
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model, Parameters
from scipy.special import erf

plt.style.use("paperdraft.mplstyle")

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

V = 20**3

def read_files_from_directory(directory):
    files = os.listdir(directory)
    files = [f for f in files if os.path.isfile(os.path.join(directory, f))]
    files.sort()
    return files

def read_file_content(file_path):
    with open(file_path, 'r') as file:
        content = file.readlines()
    values = [float(line.split()[1]) for line in content]
    return values

def extract_energy_from_filename(filename):
    match = re.search(r'E([\d.]+)sig', filename)
    if match:
        return float(match.group(1))
    return None

def fill_spectral_density_array(dir1, dir2):
    files1 = read_files_from_directory(dir1)
    files2 = read_files_from_directory(dir2)
    #assert len(files1) == len(files2), "Directories do not contain the same number of files"
    
    num_files = len(files1)
    sample_content = read_file_content(os.path.join(dir1, files1[0]))
    num_bootstraps = len(sample_content)
    spectral_density = np.zeros((num_bootstraps, 2, num_files))
    
    energy_values = []
    
    for file_index, (file1, file2) in enumerate(zip(files1, files2)):
        content1 = read_file_content(os.path.join(dir1, file1))
        content2 = read_file_content(os.path.join(dir2, file2))
        energy = extract_energy_from_filename(file1)
        energy_values.append(energy)
        
        for bootstrap_index in range(num_bootstraps):
            spectral_density[bootstrap_index][0][file_index] = content1[bootstrap_index]
            spectral_density[bootstrap_index][1][file_index] = content2[bootstrap_index]
    
    return spectral_density, energy_values

def gaussian(x, mean, sigma):
    return np.exp(-0.5 * ((x - mean) / sigma) ** 2) / (sigma * np.sqrt(np.pi/2)*(1-erf(mean/(np.sqrt(2)*sigma))))

def spectral_density_1(energy, a0, a1, E0, E1):
    sigma = 0.30
    return (a0**2 / (2 * E0) * gaussian(energy, E0, sigma) +
            a1**2 / (2 * E1) * gaussian(energy, E1, sigma))

def spectral_density_2(energy, c0, c1, E0, E1, a0, a1):
    sigma = 0.30
    return (a0 * c0 / (2 * E0) * gaussian(energy, E0, sigma) +
            a1 * c1 / (2 * E1) * gaussian(energy, E1, sigma))

def prepare_data(dir1, dir2):
    spectral_density, energy_values = fill_spectral_density_array(dir1, dir2)
    return np.array(energy_values), spectral_density

def fit_spectral_density_1(energy, spectral_density1):
    def fit_single_bootstrap(spectral_density1_sample):
        model = Model(spectral_density_1)
        params = Parameters()
        params.add('a0', value=0.000011)
        params.add('a1', value=0.000001)
        params.add('E0', value=1.0, min=0.98, max=1.02)
        params.add('E1', value=1.85, min=1.7, max=1.90)

        try:
            result = model.fit(spectral_density1_sample, params, energy=energy)
            return [result.params['a0'].value, result.params['a1'].value, result.params['E0'].value, result.params['E1'].value]
        except RuntimeError:
            return [np.nan] * 4

    fit_params = np.array([fit_single_bootstrap(sd) for sd in spectral_density1])
    
    avg_params = np.nanmean(fit_params, axis=0)
    std_params = np.nanstd(fit_params, axis=0)
    
    return avg_params, std_params, fit_params

def fit_spectral_density_2_fixed_params(energy, spectral_density2, fixed_params):
    def fit_single_bootstrap(spectral_density2_sample):
        a0, a1, E0, E1 = fixed_params
        model = Model(lambda energy, c0, c1: spectral_density_2(energy, c0, c1, E0, E1, a0, a1))
        params = Parameters()
        params.add('c0', value=0.0000010, min=0)
        params.add('c1', value=0.0000023, min=0)

        try:
            result = model.fit(spectral_density2_sample, params, energy=energy)
            return [result.params['c0'].value, result.params['c1'].value]
        except RuntimeError:
            return [np.nan] * 2

    fit_params = np.array([fit_single_bootstrap(sd) for sd in spectral_density2])
    
    avg_params = np.nanmean(fit_params, axis=0)
    std_params = np.nanstd(fit_params, axis=0)
    
    return avg_params, std_params, fit_params

def plot_with_errors(energy, avg_spectral_density1, avg_spectral_density2, fit_params_1, fit_params_2, num_bootstraps, spectral_density):
    # Increase resolution for smoother curves
    energy_fine = np.linspace(min(energy)-0.3, max(energy)+0.3, 500)
    
    # Calculate the fit curves and error bands for spectral_density_1
    fit_curves_1 = np.array([spectral_density_1(energy_fine, *params) for params in fit_params_1])
    mean_fit_curve_1 = np.mean(fit_curves_1, axis=0)
    std_fit_curve_1 = np.std(fit_curves_1, axis=0)
    
    # Calculate the fit curves and error bands for spectral_density_2
    a0, a1, E0, E1 = fit_params_1[0][:4]  # Use the first set of parameters from fit_spectral_density_1
    fit_curves_2 = np.array([spectral_density_2(energy_fine, c0, c1, E0, E1, a0, a1) for (c0, c1) in fit_params_2])
    mean_fit_curve_2 = np.mean(fit_curves_2, axis=0)
    std_fit_curve_2 = np.std(fit_curves_2, axis=0)
    
    # Error calculation for the data points
    error_spectral_density1 = np.std(spectral_density[:, 0, :], axis=0)
    error_spectral_density2 = np.std(spectral_density[:, 1, :], axis=0)

    # Calculate individual Gaussian components and their error bands for spectral_density_1
    gaussian1_component_1 = np.array([params[0]**2 / (2 * params[2]) * gaussian(energy_fine, params[2], 0.30) for params in fit_params_1])
    mean_gaussian1_component_1 = np.mean(gaussian1_component_1, axis=0)
    std_gaussian1_component_1 = np.std(gaussian1_component_1, axis=0)
    
    gaussian1_component_2 = np.array([params[1]**2 / (2 * params[3]) * gaussian(energy_fine, params[3], 0.30) for params in fit_params_1])
    mean_gaussian1_component_2 = np.mean(gaussian1_component_2, axis=0)
    std_gaussian1_component_2 = np.std(gaussian1_component_2, axis=0)

    # Calculate individual Gaussian components and their error bands for spectral_density_2
    gaussian2_component_1 = np.array([a0 * c0 / (2 * E0) * gaussian(energy_fine, E0, 0.30) for (c0, c1) in fit_params_2])
    mean_gaussian2_component_1 = np.mean(gaussian2_component_1, axis=0)
    std_gaussian2_component_1 = np.std(gaussian2_component_1, axis=0)
    
    gaussian2_component_2 = np.array([a1 * c1 / (2 * E1) * gaussian(energy_fine, E1, 0.30) for (c0, c1) in fit_params_2])
    mean_gaussian2_component_2 = np.mean(gaussian2_component_2, axis=0)
    std_gaussian2_component_2 = np.std(gaussian2_component_2, axis=0)
    
    plt.figure(figsize=(11, 5))
    
    # Plot for spectral_density_1
    plt.subplot(1, 2, 1)
    plt.errorbar(energy, avg_spectral_density1, yerr=error_spectral_density1, fmt='o', label='Data for spectral_density_1', elinewidth=1.8, markersize=6.1, markerfacecolor='none', color=CB_color_cycle[0])
    plt.plot(energy_fine, mean_fit_curve_1, color=CB_color_cycle[1], label='Mean Fit for spectral_density_1', linewidth=2.1)
    plt.fill_between(energy_fine, mean_fit_curve_1 - std_fit_curve_1, mean_fit_curve_1 + std_fit_curve_1, color=CB_color_cycle[1], alpha=0.15)
    plt.plot(energy_fine, mean_gaussian1_component_1, CB_color_cycle[5], label='Gaussian component 1')
    plt.fill_between(energy_fine, mean_gaussian1_component_1 - std_gaussian1_component_1, mean_gaussian1_component_1 + std_gaussian1_component_1, color=CB_color_cycle[5], alpha=0.15)
    plt.plot(energy_fine, mean_gaussian1_component_2, CB_color_cycle[4], label='Gaussian component 2')
    plt.fill_between(energy_fine, mean_gaussian1_component_2 - std_gaussian1_component_2, mean_gaussian1_component_2 + std_gaussian1_component_2, color=CB_color_cycle[4], alpha=0.15)
    plt.title('$N_{\\rm source} = 80$, $N_{\\rm sink} = 80$', fontsize=16)
    plt.grid(linestyle='--')
    plt.xlabel('$E/m_0$', fontsize=16)
    plt.ylabel('$\\rho_{80, 80}$', fontsize=16)
    #plt.legend()
    
    # Plot for spectral_density_2
    plt.subplot(1, 2, 2)
    plt.errorbar(energy, avg_spectral_density2, yerr=error_spectral_density2, fmt='o', label='Data for spectral_density_2', elinewidth=1.8, markersize=6.1, markerfacecolor='none', color=CB_color_cycle[0])
    plt.plot(energy_fine, mean_fit_curve_2, CB_color_cycle[2], label='Mean Fit for spectral_density_2', linewidth=2.1)
    plt.fill_between(energy_fine, mean_fit_curve_2 - std_fit_curve_2, mean_fit_curve_2 + std_fit_curve_2, color=CB_color_cycle[2], alpha=0.15)
    plt.plot(energy_fine, mean_gaussian2_component_1, CB_color_cycle[5], label='Gaussian component 1')
    plt.fill_between(energy_fine, mean_gaussian2_component_1 - std_gaussian2_component_1, mean_gaussian2_component_1 + std_gaussian2_component_1, color=CB_color_cycle[5], alpha=0.15)
    plt.plot(energy_fine, mean_gaussian2_component_2, CB_color_cycle[4], label='Gaussian component 2')
    plt.fill_between(energy_fine, mean_gaussian2_component_2 - std_gaussian2_component_2, mean_gaussian2_component_2 + std_gaussian2_component_2, color=CB_color_cycle[4], alpha=0.15)
    plt.grid(linestyle='--')
    plt.title('$N_{\\rm source} = 80$, $N_{\\rm sink} = 0$', fontsize=16)
    plt.grid(linestyle='--')
    plt.xlabel('$E/m_0$', fontsize=16)
    plt.ylabel('$\\rho_{80, 0}$', fontsize=16)
    #plt.legend()
    
    plt.tight_layout()
    # Save the plot
    plt.savefig('./spectral_density_offshell.pdf', dpi=130, bbox_inches='tight')
    #plt.show()

def main(dir1, dir2):
    mpi = 0.367
    energy, spectral_density = prepare_data(dir1, dir2)
    num_bootstraps = spectral_density.shape[0]
    energy /= mpi
    avg_spectral_density1 = np.mean(spectral_density[:, 0, :], axis=0)
    avg_spectral_density2 = np.mean(spectral_density[:, 1, :], axis=0)
    if avg_spectral_density1[3] < 0:
        spectral_density = - spectral_density
        avg_spectral_density1 = - avg_spectral_density1
        avg_spectral_density2 = - avg_spectral_density2
    
    
    avg_params_spectral_density_1, std_params_spectral_density_1, fit_params_spectral_density_1 = fit_spectral_density_1(energy, spectral_density[:, 0, :])
    
    a0, a1, E0, E1 = avg_params_spectral_density_1
    
    avg_params_spectral_density_2, std_params_spectral_density_2, fit_params_spectral_density_2 = fit_spectral_density_2_fixed_params(energy, spectral_density[:, 1, :], (a0, a1, E0, E1))
    c0, c1 = avg_params_spectral_density_2
    
    print("Parameters for spectral_density_1:")
    print(f"a0: {avg_params_spectral_density_1[0]} ± {std_params_spectral_density_1[0]}")
    print(f"a1: {avg_params_spectral_density_1[1]} ± {std_params_spectral_density_1[1]}")
    print(f"E0: {avg_params_spectral_density_1[2]} ± {std_params_spectral_density_1[2]}")
    print(f"E1: {avg_params_spectral_density_1[3]} ± {std_params_spectral_density_1[3]}")
    
    print("\nParameters for spectral_density_2 (c0, c1):")
    print(f"c0: {np.sqrt(avg_params_spectral_density_2[0] / V)} ± {std_params_spectral_density_2[0]}")
    print(f"c1: {np.sqrt(avg_params_spectral_density_2[1] / V)} ± {std_params_spectral_density_2[1]}")
    
    plot_with_errors(energy, avg_spectral_density1, avg_spectral_density2, fit_params_spectral_density_1, fit_params_spectral_density_2, num_bootstraps, spectral_density)

# Example usage
if __name__ == "__main__":
    dir1 = '/home/niccolo/work/tmp/chimera_baryons_fits/mesons/N80_N80/M3/g0g5_fund/GAUSS/g0g5_fund/Logs'  # Replace with your directory path
    dir2 = '/home/niccolo/work/tmp/chimera_baryons_fits/mesons/N80_N0/M3/g0g5_fund/GAUSS/g0g5_fund/Logs'  # Replace with your directory path
    main(dir1, dir2)
