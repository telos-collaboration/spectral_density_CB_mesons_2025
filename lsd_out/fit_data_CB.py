import argparse

import datetime
import sys
import numpy as np
import re
import matplotlib.pyplot as plt
import os
from lmfit import Parameters, Minimizer
from lsdensities.utils.rhoUtils import LogMessage
import csv
from scipy.optimize import curve_fit
import multiprocessing
import csv
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("--ensembles", nargs="+")
args = parser.parse_args()

# Plot x-limits
plot_min_lim = 0.20
plot_max_lim = 2.55

def read_csv():
    ensembles = ['M1', 'M2', 'M3', 'M4', 'M5']
    categories = ['PS', 'V', 'T', 'AV', 'AT', 'S', 'ps', 'v', 't', 'av', 'at', 's']
    Nsource_C_values_MN = {}
    Nsink_C_values_MN = {}
    am_C_values_MN = {}
    sigma1_over_mC_values_MN = {}
    sigma2_over_mC_values_MN = {}
    k_peaks = {}    # kpeaks[ensemble][channel]
    Nboot_fit = []
    with open('../metadata/metadata_spectralDensity_chimerabaryons.csv', newline='') as csvfile:
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
                Nsource_C_values_MN[ensemble].append(int(row[f'{category}_Nsource']))
                Nsink_C_values_MN[ensemble].append(int(row[f'{category}_Nsink']))
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
    categories = ['g5', 'gi', 'g0gi', 'g5gi', 'g0g5gi', 'id']
    repr = ['fund', 'as']
    ratio1 = {}
    ratio2 = {}
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ensemble = row['Ensemble']
            # Initialize lists for each ensemble if not already present
            if ensemble not in ratio1:
                ratio1[ensemble] = []
                ratio2[ensemble] = []

            for rep in repr:
                # Append data for each category to the respective lists
                for category in categories:
                    if row[f'{category}_{rep}_ratio2'] == '':
                        ratio1[ensemble].append(float(row[f'{category}_{rep}_ratio1']))
                        ratio2[ensemble].append(0.0)
                    else:
                        ratio1[ensemble].append(float(row[f'{category}_{rep}_ratio1']))
                        ratio2[ensemble].append(float(row[f'{category}_{rep}_ratio2']))
    # Create a 2D matrix with ensemble index
    matrix_2D = [
        [
            ensemble,
            ratio1[ensemble],
            ratio2[ensemble],
        ]
        for ensemble in ensembles
    ]
    return matrix_2D

def perform_fit(kernel,ensemble,rep,channel, ensemble_num, channel_num,path, file_path_input, output_name, plot_min_lim, plot_max_lim, cauchy_fit, triple_fit, four_fit, print_cov_matrix,
                plot_cov_mat, plot_corr_mat, flag_chi2, matrix_4D, k_peaks, kerneltype, nboot, fit_peaks_switch, matrix_2D):
    ####################################################################################
    '''
    # Extract directory path
    directory = os.path.dirname(output_name)
    # Check if directory exists, if not, create it
    if not os.path.exists(directory):
        os.makedirs(directory)
    '''

    # Update values if matches are found, matrix_4D[ensemble_num, object, channel_num], object = 1 --> mpi, 2 --> sigma
    mpi = matrix_4D[ensemble_num][1][channel_num]

    if kerneltype == 'GAUSS':
        sigma = matrix_4D[ensemble_num][2][channel_num]
    else:
        sigma = matrix_4D[ensemble_num][3][channel_num]
    # Define the Gaussian function
    def gaussian(x, amplitude, mean):
        return amplitude * np.exp(-((x - mean) ** 2) / (2 * sigma ** 2))

    # Define the sum of two Gaussian functions
    def double_gaussian(x, params):
        amplitude_1 = params["amplitude_1"].value
        mean_1 = params["mean_1"].value
        amplitude_2 = params["amplitude_2"].value
        mean_2 = params["mean_2"].value

        model = gaussian(x, amplitude_1, mean_1) + gaussian(x, amplitude_2, mean_2)
        return model

    def double_gaussian2(x, amplitude1, mean1, amplitude2, mean2):
        model = gaussian(x, amplitude1, mean1) + gaussian(x, amplitude2, mean2)
        return model

    def triple_gaussian2(x, amplitude1, mean1, amplitude2, mean2, amplitude3, mean3):
        model = (
                gaussian(x, amplitude1, mean1)
                + gaussian(x, amplitude2, mean2)
                + gaussian(x, amplitude3, mean3)
        )
        return model

    def four_gaussian2(
            x, amplitude1, mean1, amplitude2, mean2, amplitude3, mean3, amplitude4, mean4
    ):
        model = (
                gaussian(x, amplitude1, mean1)
                + gaussian(x, amplitude2, mean2)
                + gaussian(
            x,
            amplitude3,
            mean3,
        )
                + gaussian(x, amplitude4, mean4)
        )
        return model

    #######################################################################
    # Cauchy functions
    def cauchy(x, amplitude, mean):
        return amplitude * (sigma / ((x - mean) ** 2 + sigma ** 2))

    def double_cauchy(x, params):
        amplitude_1 = params["amplitude_1"].value
        mean_1 = params["mean_1"].value
        amplitude_2 = params["amplitude_2"].value
        mean_2 = params["mean_2"].value

        model = cauchy(x, amplitude_1, mean_1) + cauchy(x, amplitude_2, mean_2)
        return model

    def double_cauchy2(x, amplitude1, mean1, amplitude2, mean2):
        model = cauchy(x, amplitude1, mean1) + cauchy(x, amplitude2, mean2)
        return model

    def triple_cauchy2(x, amplitude1, mean1, amplitude2, mean2, amplitude3, mean3):
        model = (
                cauchy(x, amplitude1, mean1)
                + cauchy(x, amplitude2, mean2)
                + cauchy(x, amplitude3, mean3)
        )
        return model

    def four_cauchy2(
            x, amplitude1, mean1, amplitude2, mean2, amplitude3, mean3, amplitude4, mean4
    ):
        model = (
                cauchy(x, amplitude1, mean1)
                + cauchy(x, amplitude2, mean2)
                + cauchy(
            x,
            amplitude3,
            mean3,
        )
                + cauchy(x, amplitude4, mean4)
        )
        return model


    # Smearing radius in ratio with Mpi
    print(LogMessage(), "################### General Fit info #############################")
    if fit_peaks_switch == 0:
        print(LogMessage(), f"Ens: {ensemble}, Repr: {rep}, Channel: {channel}, Kernel: {kernel}, No. Peaks: {old_k_peaks}")
    elif fit_peaks_switch == 1:
        print(LogMessage(),
              f"Ens: {ensemble}, Repr: {rep}, Channel: {channel}, Kernel: {kernel}, No. Peaks: {new_k_peaks}")
    print(LogMessage(), 'sigma:', sigma)
    print(LogMessage(), 'mpi: ',mpi)
    print(LogMessage(), '##################################################################')
    # Get a list of all the file names in the directory
    file_names = os.listdir(path)
    file_names = sorted(file_names, key=lambda x: float(x.split("E")[1].split("sig")[0]))
    #print(file_names)

    # Extract the energy values from the file names
    energies = [file_name.split("E")[1].split("sig")[0] for file_name in file_names]

    # Sort the energies in ascending order
    energies.sort()
    print('energies: ', energies)
    energies = energies[:-3]
    print('energies: ', energies)
    # Define the dimensions of the matrix
    ne = len(energies)

    amplitude_vals1 = []
    mean_vals1 = []
    amplitude_vals2 = []
    mean_vals2 = []


    amplitude_vals3 = []
    mean_vals3 = []


    amplitude_vals4 = []
    mean_vals4 = []

    # Create an empty matrix
    rho_resampled = np.zeros((nboot, ne))
    #print()
    # Fill the matrix using the values from the files
    for i, energy in enumerate(energies):
        file_name = file_names[i]
        file_path = os.path.join(path, file_name)
        with open(file_path, "r") as file:
            lines = file.readlines()
            for j, line in enumerate(lines[:nboot]):
                values = line.split()
                value=float(values[1])

                if value < 0 and i < 12:
                    value = - value
                #print(value)
                rho_resampled[j, i] = float(value)

    #print(rho_resampled)
    rho_T = rho_resampled.T
    # Compute covariance matrix
    cov_matrix = np.cov(rho_T, bias=False)
    #print(cov_matrix)
    # Read the file and extract the last column of numbers
    file_path = file_path_input

    # Initialize an empty list to store the numbers from the last column
    numbers = []
    factors = np.array([])

    # Read the file and extract the last column of numbers, ignoring the first line
    with open(file_path, "r") as file:
        lines = file.readlines()
        if len(lines) > 1:  # Ensure there is at least one line in the file
            for line in lines[1:]:
                columns = line.strip().split()
                if columns:
                    # Convert the last column to a float and append it to the numbers list
                    last_column = columns[-2]
                    numbers.append(float(last_column))

    for k in range(len(np.diag(cov_matrix))):
        factors = np.append(factors, numbers[k] / np.sqrt(cov_matrix[k][k]))

    for k in range(len(np.diag(cov_matrix))):
        for j in range(len(np.diag(cov_matrix))):
            cov_matrix[j][k] *= factors[k] * factors[j]
            #cov_matrix[j][k] *= 1.0

    if print_cov_matrix is True:
        print(LogMessage(), "Evaluate covariance")
        with open(os.path.join("./covarianceMatrix_rho.txt"), "w") as output:
            for i in range(len(np.diag(cov_matrix))):
                for j in range(len(np.diag(cov_matrix))):
                    print(i, j, cov_matrix[i, j], file=output)
    '''
    print(
        LogMessage(),
        "Cond[Cov rho] = {:3.3e}".format(float(np.linalg.cond(cov_matrix))),
    )
    '''
    corrmat = np.zeros((ne, ne))
    sigmavec = cov_matrix.diagonal()

    for vi in range(ne):
        for vj in range(ne):
            corrmat[vi][vj] = cov_matrix[vi][vj] / (
                np.sqrt(sigmavec[vi]) * np.sqrt(sigmavec[vj])
            )

    if plot_cov_mat:
        plt.imshow(cov_matrix, cmap="viridis")
        plt.colorbar()
        plt.show()
    if plot_corr_mat:
        plt.imshow(corrmat, cmap="viridis")
        plt.colorbar()
        plt.show()

    np.linalg.inv(cov_matrix)

    # Extract the required columns
    x = np.array(energies, dtype=float) / mpi
    rho_central = np.zeros(ne)
    drho_central = np.zeros(ne)
    # Create a new figure with a specific size (width, height) in inches
    fig = plt.figure(figsize=(7, 4.5))  # Width: 8 inches, Height: 6 inches
    for ei in range(ne):
        rho_central[ei] = rho_T[ei].mean()
    drho_central = np.sqrt(cov_matrix.diagonal())
    # Plot the data
    plt.errorbar(
        x,
        rho_central,
        yerr=drho_central,
        fmt="o",
        color="black",
        markersize=3.0,
        label="Spectral density",
        elinewidth=1.2,
    )

    x_range = np.linspace(min(x), max(x), 1000)
    # Rough Initial guesses for parameters
    if fit_peaks_switch == 0:
        if four_fit:
            if cauchy_fit is True:
                initial_guess = [4e-7, 1.0, 6e-7, 1.5, 4e-7, 1.8, 3e-7, 2.1]
                try:
                    params, _ = curve_fit(four_cauchy2, x, rho_central, p0=initial_guess, sigma=drho_central)
                except RuntimeError as e:
                    print(f"Runtime error during fitting: {e}")
                    # Optionally, you can set default parameters or handle the error differently
                    params = initial_guess
                except Exception as e:
                    print(f"An unexpected error occurred: {e}")
                    params = initial_guess

                amp1_fit, mean1_fit, amp2_fit, mean2_fit, amp3_fit, mean3_fit, amp4_fit, mean4_fit = params
                mean1_fit = mpi
                mean2_fit = matrix_2D[ensemble_num][1][channel_num]
                mean3_fit = matrix_2D[ensemble_num][2][channel_num]
            else:
                initial_guess = [4e-7, 1.0, 6e-7, 1.5, 4e-7, 1.8, 3e-7, 2.1]
                try:
                    params, _ = curve_fit(four_gaussian2, x, rho_central, p0=initial_guess, sigma=drho_central)
                except RuntimeError as e:
                    print(f"Runtime error during fitting: {e}")
                    # Optionally, you can set default parameters or handle the error differently
                    params = initial_guess
                except Exception as e:
                    print(f"An unexpected error occurred: {e}")
                    params = initial_guess
                amp1_fit, mean1_fit, amp2_fit, mean2_fit, amp3_fit, mean3_fit, amp4_fit, mean4_fit = params
                mean1_fit = mpi
                mean2_fit = matrix_2D[ensemble_num][1][channel_num]
                mean3_fit = matrix_2D[ensemble_num][2][channel_num]
        elif (triple_fit == True and four_fit == False):
            if cauchy_fit is True:
                initial_guess = [4e-7, 1.0, 6e-7, 1.5, 4e-7, 1.9]
                try:
                    params, _ = curve_fit(triple_cauchy2, x, rho_central, p0=initial_guess, sigma=drho_central)
                except RuntimeError as e:
                    print(f"Runtime error during fitting: {e}")
                    # Optionally, you can set default parameters or handle the error differently
                    params = initial_guess
                except Exception as e:
                    print(f"An unexpected error occurred: {e}")
                    params = initial_guess

                amp1_fit, mean1_fit, amp2_fit, mean2_fit, amp3_fit, mean3_fit = params
                mean1_fit = mpi
                mean2_fit = matrix_2D[ensemble_num][1][channel_num]
                mean3_fit = matrix_2D[ensemble_num][2][channel_num]
            else:
                initial_guess = [4e-7, 1.0, 6e-7, 1.5, 4e-7, 2.1]
                try:
                    params, _ = curve_fit(triple_gaussian2, x, rho_central, p0=initial_guess, sigma=drho_central)
                except RuntimeError as e:
                    print(f"Runtime error during fitting: {e}")
                    # Optionally, you can set default parameters or handle the error differently
                    params = initial_guess
                except Exception as e:
                    print(f"An unexpected error occurred: {e}")
                    params = initial_guess

                amp1_fit, mean1_fit, amp2_fit, mean2_fit, amp3_fit, mean3_fit = params
                mean1_fit = mpi
                mean2_fit = matrix_2D[ensemble_num][1][channel_num]
                mean3_fit = matrix_2D[ensemble_num][2][channel_num]
        elif (triple_fit == False and four_fit == False):
            if cauchy_fit is True:
                initial_guess = [4e-7, 1.0, 6e-7, 1.5]
                try:
                    params, _ = curve_fit(double_cauchy2, x, rho_central, p0=initial_guess, sigma=drho_central)
                except RuntimeError as e:
                    print(f"Runtime error during fitting: {e}")
                    # Optionally, you can set default parameters or handle the error differently
                    params = initial_guess
                except Exception as e:
                    print(f"An unexpected error occurred: {e}")
                    params = initial_guess

                amp1_fit, mean1_fit, amp2_fit, mean2_fit = params
                mean1_fit = mpi
                mean2_fit = matrix_2D[ensemble_num][1][channel_num]
            else:
                initial_guess = [4e-7, 1.0, 6e-7, 1.5]
                try:
                    params, _ = curve_fit(double_gaussian2, x, rho_central, p0=initial_guess, sigma=drho_central)
                except RuntimeError as e:
                    print(f"Runtime error during fitting: {e}")
                    # Optionally, you can set default parameters or handle the error differently
                    params = initial_guess
                except Exception as e:
                    print(f"An unexpected error occurred: {e}")
                    params = initial_guess
                amp1_fit, mean1_fit, amp2_fit, mean2_fit = params
                mean1_fit = mpi
                #print('ciao: ', mean2_fit)
                mean2_fit = matrix_2D[ensemble_num][1][channel_num]
                #print('ciao2: ', mean2_fit)
    elif fit_peaks_switch == 1:
        if triple_fit == True:
            if cauchy_fit is True:
                initial_guess = [4e-7, 1.0, 6e-7, 1.5, 4e-7, 1.9]
                try:
                    params, _ = curve_fit(triple_cauchy2, x, rho_central, p0=initial_guess, sigma=drho_central)
                except RuntimeError as e:
                    print(f"Runtime error during fitting: {e}")
                    # Optionally, you can set default parameters or handle the error differently
                    params = initial_guess
                except Exception as e:
                    print(f"An unexpected error occurred: {e}")
                    params = initial_guess
                amp1_fit, mean1_fit, amp2_fit, mean2_fit, amp3_fit, mean3_fit = params
                mean1_fit = mpi
                mean2_fit = matrix_2D[ensemble_num][1][channel_num]
                mean3_fit = matrix_2D[ensemble_num][2][channel_num]
            else:
                initial_guess = [4e-7, 1.0, 6e-7, 1.5, 4e-7, 2.1]
                try:
                    params, _ = curve_fit(triple_gaussian2, x, rho_central, p0=initial_guess, sigma=drho_central)
                except RuntimeError as e:
                    print(f"Runtime error during fitting: {e}")
                    # Optionally, you can set default parameters or handle the error differently
                    params = initial_guess
                except Exception as e:
                    print(f"An unexpected error occurred: {e}")
                    params = initial_guess
                amp1_fit, mean1_fit, amp2_fit, mean2_fit, amp3_fit, mean3_fit = params
                mean1_fit = mpi
                mean2_fit = matrix_2D[ensemble_num][1][channel_num]
                mean3_fit = matrix_2D[ensemble_num][2][channel_num]
        elif (triple_fit == False and four_fit == False):
            if cauchy_fit is True:
                initial_guess = [4e-7, 1.0, 6e-7, 1.5]
                try:
                    params, _ = curve_fit(double_cauchy2, x, rho_central, p0=initial_guess, sigma=drho_central)
                except RuntimeError as e:
                    print(f"Runtime error during fitting: {e}")
                    # Optionally, you can set default parameters or handle the error differently
                    params = initial_guess
                except Exception as e:
                    print(f"An unexpected error occurred: {e}")
                    params = initial_guess
                amp1_fit, mean1_fit, amp2_fit, mean2_fit = params
                mean1_fit = mpi
                mean2_fit = matrix_2D[ensemble_num][1][channel_num]
            else:
                initial_guess = [4e-7, 1.0, 6e-7, 1.5]
                try:
                    params, _ = curve_fit(double_gaussian2, x, rho_central, p0=initial_guess, sigma=drho_central)
                except RuntimeError as e:
                    print(f"Runtime error during fitting: {e}")
                    # Optionally, you can set default parameters or handle the error differently
                    params = initial_guess
                except Exception as e:
                    print(f"An unexpected error occurred: {e}")
                    params = initial_guess
                amp1_fit, mean1_fit, amp2_fit, mean2_fit = params
                mean1_fit = mpi
                mean2_fit = matrix_2D[ensemble_num][1][channel_num]


    #print(params)
    plt.draw()

    ##################### Fitting initial parameter guesses #############################
    if fit_peaks_switch == 0:
        params = Parameters()
        params.add("amplitude_1", value=amp1_fit, min=amp1_fit - 0.4*amp1_fit, max=amp1_fit+ 0.4*amp1_fit)
        params.add("mean_1", value=1.0, min=0.99, max=1.01)
        params.add("amplitude_2", value=amp2_fit, min=amp2_fit - 0.4*amp2_fit, max=amp2_fit+ 0.4*amp2_fit)
        params.add("mean_2", value=mean2_fit, min=mean2_fit - 0.02, max=mean2_fit + 0.02)
        if triple_fit is True:
            params.add("amplitude_3", value=amp3_fit, min=amp3_fit - 0.4*amp3_fit, max=amp3_fit+ 0.4*amp3_fit)
            params.add("mean_3", value=mean3_fit, min=mean3_fit - 0.035, max=mean3_fit + 0.035)
        if four_fit is True:
            params.add("amplitude_4", value=amp4_fit, min=amp4_fit - 0.4*amp4_fit, max=amp4_fit+ 0.4*amp4_fit)
            params.add("mean_4", value=mean4_fit, min=mean4_fit - 0.4, max=mean4_fit + 0.4)
    elif fit_peaks_switch == 1:
        params = Parameters()
        params.add("amplitude_1", value=amp1_fit, min=amp1_fit - 0.4*amp1_fit, max=amp1_fit+ 0.4*amp1_fit)
        params.add("mean_1", value=1.0, min=0.99, max=1.01)
        params.add("amplitude_2", value=amp2_fit, min=amp2_fit - 0.4*amp2_fit, max=amp2_fit+ 0.4*amp2_fit)
        params.add("mean_2", value=mean2_fit, min=mean2_fit - 0.02, max=mean2_fit + 0.02)
        params.add("amplitude_3", value=1e-14, min=0.0, max=1e-10)
        params.add("mean_3", value=3.5, min=3.0, max=4.0)
        if triple_fit is True:
            params.add("amplitude_3", value=amp3_fit, min=amp3_fit - 0.4*amp3_fit, max=amp3_fit+ 0.4*amp3_fit)
            params.add("mean_3", value=mean3_fit, min=mean3_fit - 0.035, max=mean3_fit + 0.035)
            params.add("amplitude_4", value=1e-14, min=0.0, max=1e-10)
            params.add("mean_4", value=4.5, min=4.0, max=5.0)
            four_fit = True
        triple_fit = True

    #####################################################################################

    #######################################################################
    def chisq_correlated(params, x, data, cov):
        assert len(x) == len(data)

        e0 = params["mean_1"]
        e1 = params["mean_2"]
        ampl_1 = params["amplitude_1"]
        ampl_2 = params["amplitude_2"]
        if triple_fit is True:
            e2 = params["mean_3"]
            ampl_3 = params["amplitude_3"]
        if four_fit is True:
            e3 = params["mean_4"]
            ampl_4 = params["amplitude_4"]

        cov_inv = np.linalg.inv(cov)

        model = double_gaussian2(x,ampl_1, e0, ampl_2, e1)

        if cauchy_fit is True:
            model = double_cauchy2(x,ampl_1, e0, ampl_2, e1)

        if triple_fit is True:
            model = triple_gaussian2(x,ampl_1, e0, ampl_2, e1, ampl_3, e2)
            if cauchy_fit is True:
                model = triple_cauchy2(x,ampl_1, e0, ampl_2, e1, ampl_3, e2)
        if four_fit is True:
            model = four_gaussian2(x,ampl_1, e0, ampl_2, e1, ampl_3, e2, ampl_4, e3)
            if cauchy_fit is True:
                model = four_cauchy2(x,ampl_1, e0, ampl_2, e1, ampl_3, e2, ampl_4, e3)

        diff = abs(data - model)
        residual = cov_inv.dot(diff)
        return residual

    def correlated_residual(amplitude1, mean1, amplitude2, mean2, x, data, cov):
        cov_inv = np.linalg.inv(cov)
        model = double_gaussian2(x,amplitude1, mean1, amplitude2, mean2)
        if cauchy_fit is True:
            model = double_cauchy2(x,amplitude1, mean1, amplitude2, mean2)

        diff = data - model
        residual = diff * cov_inv * diff.T
        return residual

    def correlated_residual_three(
        amplitude1, mean1, amplitude2, mean2, amplitude3, mean3, x, data, cov
    ):
        cov_inv = np.linalg.inv(cov)
        model = triple_gaussian2(
            x,amplitude1, mean1, amplitude2, mean2, amplitude3, mean3
        )
        if cauchy_fit is True:
            model = triple_cauchy2(x,amplitude1, mean1, amplitude2, mean2)
        diff = model - data
        residual = diff * cov_inv * diff.T

        return residual

    def correlated_residual_four(
        amplitude1,
        mean1,
        amplitude2,
        mean2,
        amplitude3,
        mean3,
        amplitude4,
        mean4,
        x,
        data,
        cov,
    ):
        cov_inv = np.linalg.inv(cov)
        model = four_gaussian2(
            x,
            amplitude1,
            mean1,
            amplitude2,
            mean2,
            amplitude3,
            mean3,
            amplitude4,
            mean4
        )
        if cauchy_fit is True:
            model = four_cauchy2(x,amplitude1, mean1, amplitude2, mean2)
        diff = abs(model - data)
        residual = diff * cov_inv * diff.T
        return residual

    choleskyCov = np.linalg.cholesky(cov_matrix)

    for k in range(nboot):
        y = rho_resampled[k, :]

        class AltResult:
            def __init__(self, initial_guess):
                self.params = Parameters()
                self.params.add('amplitude_1', value=initial_guess[3])
                self.params.add('mean_1', value=1.0)
                self.params.add('amplitude_2', value=initial_guess[3])
                self.params.add('mean_2', value=matrix_2D[ensemble_num][1][channel_num])
                if triple_fit == True:
                    self.params.add('amplitude_3', value=initial_guess[3])
                    self.params.add('mean_3', value=matrix_2D[ensemble_num][2][channel_num])
                if four_fit == True:
                    self.params.add('amplitude_4', value=initial_guess[3])
                    self.params.add('mean_4', value=initial_guess[2])
        FITWrapper_corr = Minimizer(
            chisq_correlated, params, fcn_args=(x, y, choleskyCov)
        )
        try:
            result = FITWrapper_corr.minimize()
            fe = False
        except RuntimeError as e:
            result = AltResult(initial_guess)
            fe = True
        except Exception as e:
            result = AltResult(initial_guess)
            fe = True

        # Generate the fitted curve
        x_fit = np.linspace(plot_min_lim, plot_max_lim, 1000)

        double_gaussian(x_fit, result.params)
        if cauchy_fit is True:
            double_cauchy(x_fit, result.params)

        # Plot the fitted curve
        #        plt.plot(x_fit, y_fit, label='Fitted Curve', linewidth=1.3, color='red', alpha=0.2)

        amplitude_vals1.append(float(result.params["amplitude_1"]))
        mean_vals1.append(float(result.params["mean_1"]))
        amplitude_vals2.append(float(result.params["amplitude_2"]))
        mean_vals2.append(float(result.params["mean_2"]))

        if triple_fit is True:
            amplitude_vals3.append(float(result.params["amplitude_3"]))
            mean_vals3.append(float(result.params["mean_3"]))

        if four_fit is True:
            amplitude_vals4.append(float(result.params["amplitude_4"]))
            mean_vals4.append(float(result.params["mean_4"]))

        print(LogMessage(), "#############################")
        if fit_peaks_switch == 0:
            print(LogMessage(), f"Ens: {ensemble}, Repr: {rep}, Channel: {channel}, Kernel: {kernel}, No. Peaks: {old_k_peaks}, Bootstrap Fit number:", k, f"/ {nboot} done.")
        elif fit_peaks_switch == 1:
            print(LogMessage(),
                  f"Ens: {ensemble}, Channel: {channel}, Kernel: {kernel}, No. Peaks: {new_k_peaks}, Bootstrap Fit number:",
                  k, f"/{nboot} done.")
        '''
        print(LogMessage(), "Amplitude_1: ", float(result.params["amplitude_1"]))
        print(LogMessage(), "Mean_1: ", float(result.params["mean_1"]))
        print(LogMessage(), "Amplitude_2: ", float(result.params["amplitude_2"]))
        print(LogMessage(), "Mean_2: ", float(result.params["mean_2"]))

        if triple_fit is True:
            print(LogMessage(), "Amplitude_3: ", float(result.params["amplitude_3"]))
            print(LogMessage(), "Mean_3: ", float(result.params["mean_3"]))

        if four_fit is True:
            print(LogMessage(), "Amplitude_4: ", float(result.params["amplitude_4"]))
            print(LogMessage(), "Mean_4: ", float(result.params["mean_4"]))
        '''
        #print(LogMessage(), "#############################")
    if fe is True:
        dmean1 = 0.008*np.average(mean_vals1)
        dmean2 = 0.02*np.average(mean_vals2)
        if triple_fit is True:
            dmean3 = 0.06*np.average(mean_vals3)
        if triple_fit is True:
            dmean4 = 0.06*np.average(mean_vals4)
    else:
        dmean1 = 1.0 * np.std(mean_vals1)
        dmean2 = 2.0 * np.std(mean_vals2)
        if triple_fit is True:
            dmean3 = 0.75 * np.std(mean_vals3)
        if triple_fit is True:
            dmean4 = 0.75 * np.std(mean_vals4)
    ################## End of cycle #######################################à
    amplitude1 = np.average(amplitude_vals1)
    damplitude1 = 0.25*np.std(amplitude_vals1)
    amplitude2 = np.average(amplitude_vals2)
    damplitude2 = 0.25*np.std(amplitude_vals2)
    mean1 = np.average(mean_vals1)
    mean2 = np.average(mean_vals2)

    if triple_fit is True:
        amplitude3 = np.average(amplitude_vals3)
        damplitude3 = 0.25*np.std(amplitude_vals3)
        mean3 = np.average(mean_vals3)

    if four_fit is True:
        amplitude4 = np.average(amplitude_vals4)
        damplitude4 = 0.25*np.std(amplitude_vals4)
        mean4 = np.average(mean_vals4)


    y_gaussian_1 = [[0] * len(x_fit) for _ in range(nboot)]
    y_gaussian_2 = [[0] * len(x_fit) for _ in range(nboot)]
    y_gaussian_3 = [[0] * len(x_fit) for _ in range(nboot)]
    y_gaussian_4 = [[0] * len(x_fit) for _ in range(nboot)]
    y_gaussian_sum = [[0] * len(x_fit) for _ in range(nboot)]
    for k in range(nboot):
        # Plot the individual Gaussian components
        y_gaussian_1[k] = gaussian(x_fit, amplitude_vals1[k], mean_vals1[k])
        y_gaussian_2[k] = gaussian(x_fit, amplitude_vals2[k], mean_vals2[k])
        y_gaussian_sum[k] = y_gaussian_1[k] + y_gaussian_2[k]
        # plt.plot(x_fit, gaussian(x_fit, amplitude_vals1_new[k], mean_vals1_new[k]), label='Gaussian 1', color='orange', linewidth=1.3, alpha=0.2)
        # plt.plot(x_fit, gaussian(x_fit, amplitude_vals2_new[k], mean_vals2_new[k]), label='Gaussian 2', color='blue', linewidth=1.3, alpha=0.2)
        if cauchy_fit is True:
            y_gaussian_1[k] = cauchy(x_fit, amplitude_vals1[k], mean_vals1[k])
            y_gaussian_2[k] = cauchy(x_fit, amplitude_vals2[k], mean_vals2[k])
            y_gaussian_sum[k] = y_gaussian_1[k] + y_gaussian_2[k]

        if triple_fit is True:
            y_gaussian_3[k] = gaussian(x_fit, amplitude_vals3[k], mean_vals3[k])
            y_gaussian_sum[k] = y_gaussian_1[k] + y_gaussian_2[k] + y_gaussian_3[k]
            # plt.plot(x_fit, y_gaussian_3, label='Gaussian 3', color='gray', linewidth=1.6, alpha=0.2)
            if cauchy_fit is True:
                y_gaussian_3[k] = cauchy(x_fit, amplitude_vals3[k], mean_vals3[k])
                y_gaussian_sum[k] = y_gaussian_1[k] + y_gaussian_2[k] + y_gaussian_3[k]
        if four_fit is True:
            y_gaussian_4[k] = gaussian(x_fit, amplitude_vals4[k], mean_vals4[k])
            y_gaussian_sum[k] = (
                y_gaussian_1[k] + y_gaussian_2[k] + y_gaussian_3[k] + y_gaussian_4[k]
            )
            # plt.plot(x_fit, y_gaussian_4, label='Gaussian 4', color='brown', linewidth=1.3, alpha=0.2)
            if cauchy_fit is True:
                y_gaussian_4[k] = cauchy(x_fit, amplitude_vals4[k], mean_vals4[k])
                y_gaussian_sum[k] = (
                    y_gaussian_1[k]
                    + y_gaussian_2[k]
                    + y_gaussian_3[k]
                    + y_gaussian_4[k]
                )

    x1 = np.linspace(plot_min_lim, plot_max_lim, 1000)

    if triple_fit is True:
        if four_fit is True:
            if cauchy_fit is False:
                plt.plot(x1, gaussian(x1, amplitude4, mean4), color="mediumturquoise")
            else:
                plt.plot(x1, cauchy(x1, amplitude4, mean4), color="mediumturquoise")
            if cauchy_fit is False:
                plt.plot(
                    x1,
                    four_gaussian2(
                        x1,
                        amplitude1,
                        mean1,
                        amplitude2,
                        mean2,
                        amplitude3,
                        mean3,
                        amplitude4,
                        mean4
                    ),
                    color="gray",
                )
            else:
                plt.plot(
                    x1,
                    four_cauchy2(
                        x1,
                        amplitude1,
                        mean1,
                        amplitude2,
                        mean2,
                        amplitude3,
                        mean3,
                        amplitude4,
                        mean4,
                    ),
                    color="gray",
                )
        else:
            if cauchy_fit is False:
                plt.plot(
                    x1,
                    triple_gaussian2(
                        x1, amplitude1, mean1, amplitude2, mean2, amplitude3, mean3
                    ),
                    color="gray",
                )
            else:
                plt.plot(
                    x1,
                    triple_cauchy2(
                        x1, amplitude1, mean1, amplitude2, mean2, amplitude3, mean3
                    ),
                    color="gray",
                )
    else:
        if cauchy_fit is False:
            plt.plot(
                x1,
                double_gaussian2(x1, amplitude1, mean1, amplitude2, mean2),
                color="orange",
            )
        else:
            plt.plot(
                x1,
                double_cauchy2(x1, amplitude1, mean1, amplitude2, mean2),
                color="orange",
            )

    plt.grid(linestyle="dashed", alpha=0.6)
    # plt.legend()

    ############################# Print results ###################################################
    '''
    print("############################ Fit results ############################")
    print(
        LogMessage(),
        "E0: ",
        mean1,
        "+-",
        dmean1,
        "\t",
        "(",
        mean1 * mpi,
        "+-",
        dmean1 * mpi,
        ")",
    )
    print(
        LogMessage(),
        "E1: ",
        mean2,
        "+-",
        dmean2,
        "\t",
        "(",
        mean2 * mpi,
        "+-",
        dmean2 * mpi,
        ")",
    )

    if triple_fit is True:
        print(
            LogMessage(),
            "E2: ",
            mean3,
            "+-",
            dmean3,
            "\t",
            "(",
            mean3 * mpi,
            "+-",
            dmean3 * mpi,
            ")",
        )
    if four_fit is True:
        print(
            LogMessage(),
            "E3: ",
            mean4,
            "+-",
            dmean4,
            "\t",
            "(",
            mean4 * mpi,
            "+-",
            dmean4 * mpi,
            ")",
        )

    print(
        LogMessage(),
        "E1/E0: ",
        mean2 / mean1,
        "+-",
        np.sqrt((dmean2 / mean1) ** 2 + (dmean1 * mean2 / mean1**2) ** 2),
    )

    print(LogMessage(), "--- Fit parameters --- ")
    print(
        LogMessage(),
        "Amplitude1: ",
        amplitude1,
        "+-",
        damplitude1,
        "\t",
        "Mean1: ",
        mean1,
        "+-",
        dmean1,
    )
    print(
        LogMessage(),
        "Amplitude2: ",
        amplitude2,
        "+-",
        damplitude2,
        "\t",
        "Mean2: ",
        mean2,
        "+-",
        dmean2,
    )

    if triple_fit is True:
        print(
            LogMessage(),
            "Amplitude3: ",
            amplitude3,
            "+-",
            damplitude3,
            "\t",
            "Mean3: ",
            mean3,
            "+-",
            dmean3,
        )
    if four_fit is True:
        print(
            LogMessage(),
            "Amplitude4: ",
            amplitude4,
            "+-",
            damplitude4,
            "\t",
            "Mean4: ",
            mean4,
            "+-",
            dmean4,
        )
    '''
    # Print results to a csv
    # Print column headers
    headers = ["ensemble", "kernel", "rep", "channel", "peaks", "#aE_0", "errorE0", "#aE_1", "errorE1"]
    headers.extend(["#aE_2", "errorE2"])
    headers.extend(["#aE_3", "errorE3"])

    #print(*headers, sep=" ", file=file)
    if dmean1 < 0.001:
        dmean1 = 0.008*mean1
    if dmean2 < 0.001:
        dmean2 = 0.01*mean2
    # Print values in columns
    values = [mean1 * mpi, dmean1 * mpi, mean2 * mpi, dmean2 * mpi]

    if triple_fit:
        if dmean3 < 0.001:
            dmean3 = 0.02 * mean3
        values.extend([mean3 * mpi, dmean3 * mpi])

    if four_fit:
        if dmean4 < 0.001:
            dmean4 = 0.01 * mean4
        values.extend([mean4 * mpi, dmean4 * mpi])


    values = [
        str(value)
        for value in values
        if not isinstance(value, bool) and not isinstance(value, str)
    ]

    # Assuming 'file' is your file object opened in write mode
    with open(f'../CSVs/{ensemble}_chimerabaryons_spectral_density_spectrum.csv', 'a', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        formatted_values = [f'{float(val):.4f}' for val in values]
        if fit_peaks_switch == 0:
            csvwriter.writerow([ensemble, kernel, rep, channel, k_peaks] + formatted_values)
        elif fit_peaks_switch == 1:
            csvwriter.writerow([ensemble, kernel, rep, channel, k_peaks + 1] + formatted_values[:-2])
    result1 = gaussian(x1, amplitude1, mean1)
    transpose_y_gaussian1 = [
        [y_gaussian_1[j][i] for j in range(nboot)] for i in range(len(x1))
    ]
    if cauchy_fit is True:
        result1 = cauchy(x1, amplitude1, mean1)
        transpose_y_gaussian1 = [
            [y_gaussian_1[j][i] for j in range(nboot)] for i in range(len(x1))
        ]
    upper_band1 = [0] * len(x1)
    lower_band1 = [0] * len(x1)

    for j in range(len(x1)):
        error = 1.0 * np.std((transpose_y_gaussian1)[j])
        upper_band1[j] = result1[j] + error
        lower_band1[j] = result1[j] - error
    if cauchy_fit is False:
        plt.plot(x1, gaussian(x1, amplitude1, mean1), color="chocolate", linewidth=1.6)
        plt.fill_between(
            x1,
            lower_band1,
            upper_band1,
            color="chocolate",
            alpha=0.2,
            label="Error Bands",
        )
    else:
        plt.plot(x1, cauchy(x1, amplitude1, mean1), color="chocolate", linewidth=1.6)
        plt.fill_between(
            x1,
            lower_band1,
            upper_band1,
            color="chocolate",
            alpha=0.2,
            label="Error Bands",
        )
    result2 = gaussian(x1, amplitude2, mean2)
    transpose_y_gaussian2 = [
        [y_gaussian_2[j][i] for j in range(nboot)] for i in range(len(x1))
    ]
    if cauchy_fit is True:
        result2 = cauchy(x1, amplitude2, mean2)
        transpose_y_gaussian2 = [
            [y_gaussian_2[j][i] for j in range(nboot)] for i in range(len(x1))
        ]
    upper_band2 = [0] * len(x1)
    lower_band2 = [0] * len(x1)
    # print(transpose_y_gaussian2[500])
    for j in range(len(x1)):
        error = 1.0 * np.std((transpose_y_gaussian2)[j])
        upper_band2[j] = result2[j] + error
        lower_band2[j] = result2[j] - error

    if cauchy_fit is False:
        plt.plot(x1, gaussian(x1, amplitude2, mean2), color="olive", linewidth=1.6)
        plt.fill_between(
            x1, lower_band2, upper_band2, color="olive", alpha=0.25, label="Error Bands"
        )
    else:
        plt.plot(x1, cauchy(x1, amplitude2, mean2), color="olive", linewidth=1.6)
        plt.fill_between(
            x1, lower_band2, upper_band2, color="olive", alpha=0.25, label="Error Bands"
        )
    if four_fit is True:
        if cauchy_fit is False:
            plt.plot(x1, gaussian(x1, amplitude3, mean3), color="orange", linewidth=1.6)

            result3 = gaussian(x1, amplitude3, mean3)
            transpose_y_gaussian3 = [
                [y_gaussian_3[j][i] for j in range(nboot)] for i in range(len(x1))
            ]
            upper_band3 = [0] * len(x1)
            lower_band3 = [0] * len(x1)
            for j in range(len(x1)):
                error = 1.0 * np.std((transpose_y_gaussian3)[j])
                upper_band3[j] = result3[j] + error
                lower_band3[j] = result3[j] - error

            plt.fill_between(
                x1,
                lower_band3,
                upper_band3,
                color="orange",
                alpha=0.25,
                label="Error Bands",
            )

            result_sum = (
                gaussian(x1, amplitude2, mean2)
                + gaussian(x1, amplitude1, mean1)
                + gaussian(x1, amplitude3, mean3)
                + gaussian(x1, amplitude4, mean4)
            )

            transpose_y_gaussian_sum = [
                [y_gaussian_sum[j][i] for j in range(nboot)] for i in range(len(x1))
            ]
            upper_band5 = [0] * len(x1)
            lower_band5 = [0] * len(x1)
            # print(transpose_y_gaussian1[999])
            for j in range(len(x1)):
                error = 1.0 * np.std((transpose_y_gaussian_sum)[j])
                upper_band5[j] = result_sum[j] + error
                lower_band5[j] = result_sum[j] - error

            plt.plot(x1, gaussian(x1, amplitude4, mean4), color="pink", linewidth=1.6)

            result4 = gaussian(x1, amplitude4, mean4)
            transpose_y_gaussian4 = [
                [y_gaussian_4[j][i] for j in range(nboot)] for i in range(len(x1))
            ]
            upper_band4 = [0] * len(x1)
            lower_band4 = [0] * len(x1)
            for j in range(len(x1)):
                error = 1.0 * np.std((transpose_y_gaussian4)[j])
                upper_band4[j] = result4[j] + error
                lower_band4[j] = result4[j] - error

            plt.fill_between(
                x1,
                lower_band4,
                upper_band4,
                color="pink",
                alpha=0.25,
                label="Error Bands",
            )

        else:
            plt.plot(x1, cauchy(x1, amplitude3, mean3), color="orange", linewidth=1.6)

            result3 = cauchy(x1, amplitude3, mean3)
            transpose_y_gaussian3 = [
                [y_gaussian_3[j][i] for j in range(nboot)] for i in range(len(x1))
            ]
            upper_band3 = [0] * len(x1)
            lower_band3 = [0] * len(x1)
            for j in range(len(x1)):
                error = 1.0 * np.std((transpose_y_gaussian3)[j])
                upper_band3[j] = result3[j] + error
                lower_band3[j] = result3[j] - error

            plt.fill_between(
                x1,
                lower_band3,
                upper_band3,
                color="orange",
                alpha=0.25,
                label="Error Bands",
            )

            result_sum = (
                cauchy(x1, amplitude2, mean2)
                + cauchy(x1, amplitude1, mean1)
                + cauchy(x1, amplitude3, mean3)
                + cauchy(x1, amplitude4, mean4)
            )

            transpose_y_gaussian_sum = [
                [y_gaussian_sum[j][i] for j in range(nboot)] for i in range(len(x1))
            ]
            upper_band5 = [0] * len(x1)
            lower_band5 = [0] * len(x1)
            # print(transpose_y_gaussian1[999])
            for j in range(len(x1)):
                error = 1.0 * np.std((transpose_y_gaussian_sum)[j])
                upper_band5[j] = result_sum[j] + error
                lower_band5[j] = result_sum[j] - error

            plt.plot(x1, cauchy(x1, amplitude4, mean4), color="pink", linewidth=1.6)

            result4 = cauchy(x1, amplitude4, mean4)
            transpose_y_gaussian4 = [
                [y_gaussian_4[j][i] for j in range(nboot)] for i in range(len(x1))
            ]
            upper_band4 = [0] * len(x1)
            lower_band4 = [0] * len(x1)
            for j in range(len(x1)):
                error = 1.0 * np.std((transpose_y_gaussian4)[j])
                upper_band4[j] = result4[j] + error
                lower_band4[j] = result4[j] - error

            plt.fill_between(
                x1,
                lower_band4,
                upper_band4,
                color="pink",
                alpha=0.25,
                label="Error Bands",
            )

        plt.fill_between(
            x1, lower_band5, upper_band5, color="gray", alpha=0.25, label="Error Bands"
        )
    elif triple_fit is True:
        if cauchy_fit is False:
            plt.plot(x1, gaussian(x1, amplitude3, mean3), color="orange", linewidth=1.6)

            result_sum = (
                gaussian(x1, amplitude2, mean2)
                + gaussian(x1, amplitude1, mean1)
                + gaussian(x1, amplitude3, mean3)
            )
            result3 = gaussian(x1, amplitude3, mean3)
            transpose_y_gaussian3 = [
                [y_gaussian_3[j][i] for j in range(nboot)] for i in range(len(x1))
            ]
            upper_band3 = [0] * len(x1)
            lower_band3 = [0] * len(x1)
            for j in range(len(x1)):
                error = 1.0 * np.std((transpose_y_gaussian3)[j])
                upper_band3[j] = result3[j] + error
                lower_band3[j] = result3[j] - error

            plt.fill_between(
                x1,
                lower_band3,
                upper_band3,
                color="orange",
                alpha=0.25,
                label="Error Bands",
            )

        else:
            plt.plot(x1, cauchy(x1, amplitude3, mean3), color="orange", linewidth=1.6)

            result_sum = (
                cauchy(x1, amplitude2, mean2)
                + cauchy(x1, amplitude1, mean1)
                + cauchy(x1, amplitude3, mean3)
            )

            result3 = cauchy(x1, amplitude3, mean3)
            transpose_y_gaussian3 = [
                [y_gaussian_3[j][i] for j in range(nboot)] for i in range(len(x1))
            ]
            upper_band3 = [0] * len(x1)
            lower_band3 = [0] * len(x1)
            for j in range(len(x1)):
                error = 1.0 * np.std((transpose_y_gaussian3)[j])
                upper_band3[j] = result3[j] + error
                lower_band3[j] = result3[j] - error

            plt.fill_between(
                x1,
                lower_band3,
                upper_band3,
                color="orange",
                alpha=0.25,
                label="Error Bands",
            )

        transpose_y_gaussian_sum = [
            [y_gaussian_sum[j][i] for j in range(nboot)] for i in range(len(x1))
        ]
        upper_band4 = [0] * len(x1)
        lower_band4 = [0] * len(x1)
        # print(transpose_y_gaussian1[999])
        for j in range(len(x1)):
            error = 1.0 * np.std((transpose_y_gaussian_sum)[j])
            upper_band4[j] = result_sum[j] + error
            lower_band4[j] = result_sum[j] - error

        plt.fill_between(
            x1, lower_band4, upper_band4, color="gray", alpha=0.25, label="Error Bands"
        )
    else:
        if cauchy_fit is False:
            plt.plot(
                x1,
                double_gaussian2(x1, amplitude1, mean1, amplitude2, mean2),
                color="orange",
                linewidth=1.8,
            )
            result_sum = gaussian(x1, amplitude2, mean2) + gaussian(
                x1, amplitude1, mean1
            )
        else:
            plt.plot(
                x1,
                double_cauchy2(x1, amplitude1, mean1, amplitude2, mean2),
                color="orange",
                linewidth=1.8,
            )
            result_sum = cauchy(x1, amplitude2, mean2) + cauchy(x1, amplitude1, mean1)

        transpose_y_gaussian_sum = [
            [y_gaussian_sum[j][i] for j in range(nboot)] for i in range(len(x1))
        ]
        upper_band3 = [0] * len(x1)
        lower_band3 = [0] * len(x1)
        # print(transpose_y_gaussian1[999])
        for j in range(len(x1)):
            error = 1.0 * np.std((transpose_y_gaussian_sum)[j])
            upper_band3[j] = result_sum[j] + error
            lower_band3[j] = result_sum[j] - error

        plt.fill_between(
            x1,
            lower_band3,
            upper_band3,
            color="orange",
            alpha=0.25,
            label="Error Bands",
        )

    if flag_chi2:
        if triple_fit is True:
            if four_fit is True:
                if len(x) > 8:
                    chi_square_red = correlated_residual_four(
                        amplitude1,
                        mean1,
                        amplitude2,
                        mean2,
                        amplitude3,
                        mean3,
                        amplitude4,
                        mean4,
                        x,
                        rho_central,
                        cov_matrix,
                    ) / (len(x) - 8)
                else:
                    print("Cannot compute Chi square!")
                    flag_chi2 = False
            else:
                if len(x) > 6:
                    chi_square_red = correlated_residual_three(
                        amplitude1,
                        mean1,
                        amplitude2,
                        mean2,
                        amplitude3,
                        mean3,
                        x,
                        rho_central,
                        cov_matrix,
                    )[0, 0] / (len(x) - 6)
                else:
                    print("Cannot compute Chi square!")
                    flag_chi2 = False
        else:
            if len(x) > 4:
                chi_square_red = correlated_residual(
                    amplitude1, mean1, amplitude2, mean2, x, rho_central, cov_matrix
                )[0, 0] / (len(x) - 4)
                print("len(x): ", len(x))
            else:
                print("Cannot compute Chi square!")
                flag_chi2 = False
        print(LogMessage(), " Reduced Chi Square (with correlation): ", chi_square_red)

    # Plot the data
    plt.errorbar(
        x,
        rho_central,
        yerr=drho_central,
        fmt="o",
        color="black",
        markersize=3.0,
        label="Spectral density",
        elinewidth=1.2,
    )
    '''
    # Save the figure with the specified size
    plt.savefig(output_name, format="pdf", dpi=300)
    '''
    # Display the plot
    #plt.show()
    plt.close(fig)
    return None

########################### Preferences ################################
# If you want to fit with Cauchy (False == Gaussians)
cauchy_fit = True
# If both false, it's two Gaussians/Cauchy fit
triple_fit = False
four_fit = False
print_cov_matrix = False
plot_cov_mat = False
plot_corr_mat = False
flag_chi2 = False  # To compute and show the correlated chi-square
if four_fit is True:
    triple_fit = True

matrix_4D, k_peaks, Nboot_fit  = read_csv()
file_path_MD = '../metadata/ratioguesses_chimerabaryons_spectrum.csv'
matrix_2D = read_csv2(file_path_MD)
ensembles = ['M1', 'M2', 'M3', 'M4', 'M5']
#ensembles = ['M1']
mesonic_channels = ['Chimera_OC_even', 'Chimera_OV12_even', 'Chimera_OV32_even', 'Chimera_OC_odd',  'Chimera_OV12_odd', 'Chimera_OV32_odd']
#mesonic_channels = ['id']
reps = ['as']
#reps = ['as']
kerneltype = ['GAUSS', 'CAUCHY']
#kerneltype = ['GAUSS']
#ensemble_num = 4
#channel_num = 5
Nsource = 80
Nsink = 0

headers = ["Label", "kernel", "rep", "channel", "peaks", "aE_0", "errorE0", "aE_1", "errorE1"]
headers.extend(["aE_2", "errorE2"])
headers.extend(["aE_3", "errorE3"])


# TODO: match names with spec_dens code outputs in our inputs
# TODO: sp_dens_code.py --> structure of 'input_fit/'

for ensemble in args.ensembles:
    with open(f'../CSVs/{ensemble}_chimerabaryons_spectral_density_spectrum.csv', 'a', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(headers)
    ensemble_num = ensembles.index(ensemble)
    for rep in reps:
        for k, channel in enumerate(mesonic_channels):
            for kernel in kerneltype:
                for fit_peaks_switch in range(2):
                    new_k_peaks = k_peaks[ensemble][k] + 1
                    old_k_peaks = k_peaks[ensemble][k]
                    channel_num = k
                    if rep == 'as':
                        channel_num += 6

                    if kernel == 'GAUSS':
                        cauchy_fit = False
                    elif kernel == 'CAUCHY':
                        cauchy_fit = True
                    if k_peaks[ensemble][channel_num] == 2:  # kpeaks[ensemble][channel]
                        triple_fit = False
                        four_fit = False
                    elif k_peaks[ensemble][channel_num] == 3:
                        triple_fit = True
                        four_fit = False
                    elif k_peaks[ensemble][channel_num] == 4:
                        triple_fit = True
                        four_fit = True

                    path = f"../input_fit/{ensemble}/{channel}_Nsource{Nsource}_Nsink{Nsink}/{kernel}/{channel}_Nsource{Nsource}_Nsink{Nsink}/Logs/"
                    file_path_input = f"../input_fit/{ensemble}/{channel}_Nsource{Nsource}_Nsink{Nsink}/{kernel}/fit_results.txt"

                    if fit_peaks_switch == 0:
                        output_name = f"./fitresults/fit_results_{ensemble}_{channel}_{kernel}_kpeaks{old_k_peaks}.pdf"
                    elif fit_peaks_switch == 1:
                        output_name = f"./fitresults/fit_results_{ensemble}_{channel}_{kernel}_kpeaks{new_k_peaks}.pdf"

                    perform_fit(kernel,ensemble,rep,channel,ensemble_num, channel_num, path, file_path_input, output_name, plot_min_lim, plot_max_lim, cauchy_fit, triple_fit, four_fit, print_cov_matrix,
                                    plot_cov_mat, plot_corr_mat, flag_chi2, matrix_4D, k_peaks[ensemble][channel_num], kernel, Nboot_fit[ensemble_num], fit_peaks_switch, matrix_2D)




# Avoid needing to work out the full tangle of output files,
# while still allowing a workflow dependency on completing this rule
with open("fit_data_CB_complete", "w") as completion_tag_file:
    print(
        f"CB fitting complete at {datetime.datetime.now().astimezone()}",
        file=completion_tag_file,
    )
