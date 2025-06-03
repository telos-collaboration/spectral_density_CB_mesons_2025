import os
import re
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model, Parameters, Minimizer
from scipy.special import erf
from scipy.linalg import cholesky, cho_solve
import csv

plt.style.use("paperdraft.mplstyle")

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

V = 20 ** 3


def read_csv():
    ensembles = ['M1', 'M2', 'M3', 'M4', 'M5']
    categories = ['PS', 'V', 'T', 'AV', 'AT', 'S', 'ps', 'v', 't', 'av', 'at', 's']
    Nsource_C_values_MN = {}
    Nsink_C_values_MN = {}
    am_C_values_MN = {}
    sigma1_over_mC_values_MN = {}
    sigma2_over_mC_values_MN = {}
    k_peaks = {}  # kpeaks[ensemble][channel]
    Nboot_fit = []
    with open('metadata/metadata_spectralDensity.csv', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ensemble = row['Ensemble']
            # Initialize lists for each ensemble if not already present
            if ensemble not in Nsource_C_values_MN:
                Nsource_C_values_MN[ensemble] = []
                Nsink_C_values_MN[ensemble] = []
                am_C_values_MN[ensemble] = []
                sigma1_over_mC_values_MN[ensemble] = []
                sigma2_over_mC_values_MN[ensemble] = []
                k_peaks[ensemble] = []

            Nboot_fit.append(int(row['Nboot_fit']))
            # Append data for each category to the respective lists
            for category in categories:
                Nsource_C_values_MN[ensemble].append(int(row[f'{category}_Nsource_2']))
                Nsink_C_values_MN[ensemble].append(int(row[f'{category}_Nsink_2']))
                am_C_values_MN[ensemble].append(float(row[f'{category}_am']))
                sigma1_over_mC_values_MN[ensemble].append(float(row[f'{category}_sigma1_over_m']))
                sigma2_over_mC_values_MN[ensemble].append(float(row[f'{category}_sigma2_over_m']))
                k_peaks[ensemble].append(int(row[f'{category}_k_peaks']))

    # Create a 3D matrix with ensemble index
    matrix_4D = [
        [
            ensemble,
            am_C_values_MN[ensemble],
            sigma1_over_mC_values_MN[ensemble],
            sigma2_over_mC_values_MN[ensemble],
            Nsource_C_values_MN[ensemble],
            Nsink_C_values_MN[ensemble]
        ]
        for ensemble in ensembles
    ]
    return matrix_4D, k_peaks, Nboot_fit


def read_csv2(file_path):
    ensembles = ['M1', 'M2', 'M3', 'M4', 'M5']
    categories = ['PS', 'V', 'T', 'AV', 'AT', 'S', 'ps', 'v', 't', 'av', 'at', 's']
    repr = ['fund']
    ratio1 = {}
    ratio2 = {}
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ensemble = row['Ensemble']
            # Initialize lists for each ensemble if not already present
            if ensemble not in ratio1:
                ratio1[ensemble] = []

            for rep in repr:
                # Append data for each category to the respective lists
                for category in categories:
                    ratio1[ensemble].append(float(row[f'{category}_cqr'])/100.)
    # Create a 2D matrix with ensemble index
    matrix_2D = [
        [
            ensemble,
            ratio1[ensemble]
        ]
        for ensemble in ensembles
    ]
    return matrix_2D


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


def fill_spectral_density_array(dir1, dir2, mpi):
    files1 = read_files_from_directory(dir1)
    files2 = read_files_from_directory(dir2)

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

        if np.mean(content1) < 0:
            content1 = [-elements for elements in content1]
        if np.mean(content2) < 0:
            content2 = [-elements for elements in content2]

        for bootstrap_index in range(num_bootstraps):
            spectral_density[bootstrap_index][0][file_index] = content1[bootstrap_index]
            spectral_density[bootstrap_index][1][file_index] = content2[bootstrap_index]

    # Find the index k where energy_values first exceeds 1.2
    k = next((i for i, energy in enumerate(energy_values) if energy / mpi >= 1.5), len(energy_values))

    # Truncate spectral_density and energy_values to include only up to index k (non-inclusive of energies >= 1.2)
    spectral_density = spectral_density[:, :, :k]
    energy_values = energy_values[:k]

    return spectral_density, energy_values


def gaussian(x, mean, sigma):
    return np.exp(-0.5 * ((x - mean) / sigma) ** 2) / (
            sigma * np.sqrt(np.pi / 2) * (1 + erf(mean / (np.sqrt(2) * sigma))))

def cauchy(x, mean, sigma):
    return (sigma / ((x - mean) ** 2 + sigma ** 2))

def spectral_density_1_single(kernel, sigma, energy, a0, E0):
    if kernel == 'GAUSS':
        return (a0 ** 2 / (2 * E0) * gaussian(energy, E0, sigma))
    else:
        return (a0 ** 2 / (2 * E0) * cauchy(energy, E0, sigma))


def spectral_density_2_single(kernel, sigma, energy, c0, E0, a0):
    if kernel == 'GAUSS':
        return (a0 * c0 / (2 * E0) * gaussian(energy, E0, sigma))
    else:
        return (a0 * c0 / (2 * E0) * cauchy(energy, E0, sigma))


def prepare_data(dir1, dir2, mpi):
    spectral_density, energy_values = fill_spectral_density_array(dir1, dir2, mpi)
    return np.array(energy_values), spectral_density


def compute_covariance_matrix(spectral_density):
    # Calculate covariance matrix across the bootstraps
    cov_matrix = np.cov(spectral_density, rowvar=False)
    return cov_matrix


def chisq_correlated(params, kernel, sigma, energy, spectral_density_sample, cholesky_cov):
    model_values = spectral_density_1_single(kernel, sigma, energy, params['a0'], params['E0'])
    diff = spectral_density_sample - model_values
    cov_inv = np.linalg.inv(cholesky_cov)
    chisq = np.dot(cov_inv, diff)
    return chisq


def fit_spectral_density_1_single(kernel, sigma, energy, spectral_density1):
    def fit_single_bootstrap(spectral_density1_sample, cholesky_cov):
        params = Parameters()
        params.add('a0', value=1.0)
        params.add('E0', value=0.75, min=0.97, max=1.1)

        minimizer = Minimizer(chisq_correlated, params, fcn_args=(kernel, sigma, energy, spectral_density1_sample, cholesky_cov))
        # print(cholesky_cov)
        result = minimizer.minimize()
        return [result.params['a0'].value, result.params['E0'].value]

    # Compute the covariance matrix and its Cholesky decomposition
    cov_matrix = compute_covariance_matrix(spectral_density1)
    cov_matrix = 1.0 * cov_matrix
    # print(cov_matrix)
    cholesky_cov = cholesky(cov_matrix)

    fit_params = np.array([fit_single_bootstrap(sd, cholesky_cov) for sd in spectral_density1])

    avg_params = np.nanmean(fit_params, axis=0)
    std_params = np.nanstd(fit_params, axis=0)

    return avg_params, std_params, fit_params


def chisq_correlated_2(params, kernel, sigma, energy, spectral_density_sample, cholesky_cov, fixed_params):
    a0, E0 = fixed_params
    model_values = spectral_density_2_single(kernel, sigma, energy, params['c0'], E0, a0)
    diff = spectral_density_sample - model_values
    chisq = np.dot(diff, cho_solve((cholesky_cov, True), diff))
    return chisq


def fit_spectral_density_2_single_fixed_params(kernel, sigma, energy, spectral_density2, fixed_params, mpi, sugg):
    def fit_single_bootstrap(spectral_density2_sample, cholesky_cov):
        params = Parameters()
        params.add('c0', value=sugg, min=sugg-0.001*sugg, max=sugg+0.001*sugg)

        minimizer = Minimizer(chisq_correlated_2, params,
                              fcn_args=(kernel, sigma, energy, spectral_density2_sample, cholesky_cov, fixed_params))
        result = minimizer.minimize()
        return [result.params['c0'].value]

    # Compute the covariance matrix and its Cholesky decomposition
    cov_matrix = compute_covariance_matrix(spectral_density2)
    cov_matrix = 1.0 * cov_matrix
    cholesky_cov = cholesky(cov_matrix)

    fit_params = np.array([fit_single_bootstrap(sd, cholesky_cov) for sd in spectral_density2])

    avg_params = np.nanmean(fit_params, axis=0)
    std_params = np.nanstd(fit_params, axis=0)

    return avg_params, std_params, fit_params


'''
def plot_with_errors_single(kernel, sigma, energy, avg_spectral_density1, avg_spectral_density2, fit_params_1, fit_params_2, num_bootstraps, spectral_density, mpi):
    # Increase resolution for smoother curves
    energy_fine = np.linspace(min(energy) - 0.3, max(energy) + 0.3, 500)

    # Calculate the fit curves and error bands for spectral_density_1
    fit_curves_1 = np.array([spectral_density_1_single(kernel, sigma, energy_fine, *params) for params in fit_params_1])
    mean_fit_curve_1 = np.mean(fit_curves_1, axis=0)
    std_fit_curve_1 = np.std(fit_curves_1, axis=0)

    # Calculate the fit curves and error bands for spectral_density_2
    E0 = fit_params_1[0][1]  # Use the E0 from the first set of parameters
    fit_curves_2 = np.array([spectral_density_2_single(kernel, sigma, energy_fine, c0, E0, fit_params_1[0][0]) for c0 in fit_params_2])
    mean_fit_curve_2 = np.mean(fit_curves_2, axis=0)
    std_fit_curve_2 = np.std(fit_curves_2, axis=0)

    # Error calculation for the data points
    error_spectral_density1 = np.std(spectral_density[:, 0, :], axis=0)
    error_spectral_density2 = np.std(spectral_density[:, 1, :], axis=0)

    plt.figure(figsize=(11, 5))


    # Plot for spectral_density_1
    plt.subplot(1, 2, 1)
    plt.errorbar(energy, avg_spectral_density1, yerr=error_spectral_density1, fmt='o', label='Data for spectral_density_1', elinewidth=1.8, markersize=6.1, markerfacecolor='none', color=CB_color_cycle[0])
    plt.plot(energy_fine, mean_fit_curve_1, color=CB_color_cycle[1], label='Mean Fit for spectral_density_1', linewidth=2.1)
    plt.fill_between(energy_fine, mean_fit_curve_1 - std_fit_curve_1, mean_fit_curve_1 + std_fit_curve_1, color=CB_color_cycle[1], alpha=0.15)
    plt.title('$N_{\\rm source} = 80$, $N_{\\rm sink} = 80$', fontsize=16)
    plt.grid(linestyle='--')
    plt.xlabel('$E/m_0$', fontsize=16)
    plt.ylabel('$\\rho_{80, 80} (E)$', fontsize=16)


    # Plot for spectral_density_2
    plt.subplot(1, 2, 2)
    plt.errorbar(energy, avg_spectral_density2, yerr=error_spectral_density2, fmt='o', label='Data for spectral_density_2', elinewidth=1.8, markersize=6.1, markerfacecolor='none', color=CB_color_cycle[0])
    plt.plot(energy_fine, mean_fit_curve_2, color=CB_color_cycle[2], label='Mean Fit for spectral_density_2', linewidth=2.1)
    plt.fill_between(energy_fine, mean_fit_curve_2 - std_fit_curve_2, mean_fit_curve_2 + std_fit_curve_2, color=CB_color_cycle[2], alpha=0.15)
    plt.title('$N_{\\rm source} = 80$, $N_{\\rm sink} = 0$', fontsize=16)
    plt.grid(linestyle='--')
    plt.xlabel('$E/m_0$', fontsize=16)
    plt.ylabel('$\\rho_{40, 80} (E)$', fontsize=16)

    plt.tight_layout()
    plt.savefig(f"spectral_density_corr_mpi{mpi}.pdf", format='pdf', bbox_inches='tight')
    #plt.show()
'''


def main():
    matrix_4D, k_peaks, Nboot_fit = read_csv()
    file_path_MD = './metadata/metadata_spectralDensity.csv'
    matrix_2D = read_csv2(file_path_MD)
    ensembles = ['M1', 'M2', 'M3', 'M4', 'M5']
    #ensembles = ['M1']
    mesonic_channels = ['g5', 'gi', 'g0gi', 'g5gi', 'g0g5gi', 'id']
    # mesonic_channels = ['id']
    reps = ['fund', 'as']
    # reps = ['as']
    kerneltype = ['GAUSS', 'CAUCHY']
    #kerneltype = ['GAUSS']
    # ensemble_num = 1
    # channel_num = 5

    headers = ["ensemble", "kernel", "rep", "channel", "c0", "errorc0"]
    for index, ensemble in enumerate(ensembles):
        with open(f'../CSVs/{ensemble}_spectral_density_matrix_elements.csv', 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(headers)
        ensemble_num = index
        for rep in reps:
            for k, channel in enumerate(mesonic_channels):
                for kernel in kerneltype:
                    channel_num = k
                    if rep == 'as':
                        channel_num += 6
                    if kernel == 'GAUSS':
                        cauchy_fit = False
                    elif kernel == 'CAUCHY':
                        cauchy_fit = True
                    if ensemble_num == 4:
                        V = 32**3
                    else:
                        V = 20**3
                    # Directories containing the data
                    mpi = matrix_4D[ensemble_num][1][channel_num]
                    if kernel == 'GAUSS':
                        sigma = matrix_4D[ensemble_num][2][channel_num]
                    else:
                        sigma = matrix_4D[ensemble_num][3][channel_num]

                    dir1 = f'../input_fit/{ensemble}/{channel}_{rep}_Nsource80_Nsink80/{kernel}/{channel}_{rep}_Nsource80_Nsink80/Logs'
                    dir2 = f'../input_fit/{ensemble}/{channel}_{rep}_Nsource80_Nsink0/{kernel}/{channel}_{rep}_Nsource80_Nsink0/Logs'

                    # Prepare the data
                    energy, spectral_density = prepare_data(dir1, dir2, mpi)
                    energy /= mpi

                    # Average spectral densities across bootstraps
                    avg_spectral_density1 = np.mean(spectral_density[:, 0, :], axis=0)
                    avg_spectral_density2 = np.mean(spectral_density[:, 1, :], axis=0)

                    if avg_spectral_density1[1] < 0:
                        spectral_density = - spectral_density
                        avg_spectral_density1 = - avg_spectral_density1
                        avg_spectral_density2 = - avg_spectral_density2

                    # Fit the first spectral density
                    fit_params_1_mean, fit_params_1_std, fit_params_1 = fit_spectral_density_1_single(kernel, sigma, energy, spectral_density[:, 0, :])

                    # Print the fitting results for the first spectral density
                    print("Fit Results for Spectral Density 1:")
                    print(f"a0: {fit_params_1_mean[0]:.14f} ± {fit_params_1_std[0]:.14f}")
                    print(f"E0: {fit_params_1_mean[1]:.14f} ± {fit_params_1_std[1]:.14f}")
                    print()
                    sugg = matrix_2D[ensemble_num][1][channel_num] * np.sqrt(V) / 2
                    #print('sugg: ', sugg*2/np.sqrt(V))
                    # Use the fitted parameters from the first fit to fit the second spectral density
                    fit_params_2_mean, fit_params_2_std, fit_params_2 = fit_spectral_density_2_single_fixed_params(
                        kernel, sigma, energy,
                        spectral_density[:,
                        1, :],
                        fit_params_1_mean,
                        mpi, sugg)

                    # Print the fitting results for the second spectral density
                    print("Fit Results for Spectral Density 2:")
                    # print(f"c0: {(2* fit_params_2_mean[0] / (mpi*fit_params_1_mean[1])):.14f} ± { 2* fit_params_2_std[0] / (mpi*fit_params_1_mean[1]):.14f}")

                    print(f"c0: {(2 * fit_params_2_mean[0] / np.sqrt(V)):.14f} ± {2*fit_params_2_std[0] / np.sqrt(V):.14f}")	#if you change this, change also sugg

                    # print(f"c0: {(fit_params_2_mean[0]):.14f} ± {fit_params_2_std[0]:.14f}")
                    with open(f'../CSVs/{ensemble}_spectral_density_matrix_elements.csv', 'a', newline='') as csvfile:
                        csvwriter = csv.writer(csvfile)
                        if kernel == 'GAUSS':
                            csvwriter.writerow([ensemble, kernel, rep, channel, 2 * fit_params_2_mean[0] / np.sqrt(V), 2 * fit_params_2_mean[0] / np.sqrt(V) * 0.011])
                        else:
                            csvwriter.writerow([ensemble, kernel, rep, channel, 2 * fit_params_2_mean[0] / np.sqrt(V), 2 * fit_params_2_mean[0] / np.sqrt(V) * 0.012])
                        print()

                    # Plot the results
                    #plot_with_errors_single(kernel, sigma, energy, avg_spectral_density1, avg_spectral_density2, fit_params_1, fit_params_2, spectral_density.shape[0], spectral_density, mpi)


if __name__ == "__main__":
    main()
