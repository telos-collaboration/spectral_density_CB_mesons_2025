import pandas as pd
import numpy as np
import json
import math


def read_json(file_path):
    """Reads a JSON file and returns the parsed data."""
    with open(file_path, "r") as file:
        data = json.load(file)
    return data


def bootstrap_errors(data, num_resamples=1000):
    """Calculates the bootstrap error from the data."""
    resamples = np.random.choice(data, size=(num_resamples, len(data)), replace=True)
    means = np.mean(resamples, axis=1)
    mean_of_means = np.mean(means)
    std_error = np.std(means)
    return mean_of_means, std_error


def process_gevp_En_mass_samples(file_path, channel, rep, n):
    """Processes the 'gevp_f_at_E0_mass_samples' from the JSON file and calculates the average and bootstrap errors."""
    # Read JSON file
    data = read_json(file_path)
    print(f"gevp_{rep}_{channel}_E{n}_mass_samples")
    # Extract 'gevp_f_at_E0_mass_samples'
    samples = np.array(data[f"gevp_{rep}_{channel}_E{n}_mass_samples"])

    # Calculate the average and bootstrap error
    average, error = bootstrap_errors(samples)

    return average, error


def add_error(channel_E0, err):
    if channel_E0 != 0 and channel_E0 != "NaN" and channel_E0 != "-":
        # print(err)
        # print('channel_E0_2: ', channel_E0)
        # print('err_channel_E0_2: ' , err)
        # Convert the error to a string with significant digits
        err_str = f"{err:.2g}".replace(
            ".", ""
        )  # Convert error to string with 2 significant digits and remove decimal
        # print(err_str)
        # Count leading zeros
        leading_zeros = len(err_str) - len(err_str.lstrip("0"))

        # Remove leading zeros
        err_str = err_str.lstrip("0")

        # print(err_str)
        if len(str(err_str)) == 1:
            err_str = str(int(err_str) * 30)
        # Determine the number of significant digits in the error
        err_significant_digits = len(err_str) + leading_zeros

        # Format the channel_E0 value with the correct number of significant digits
        # Calculate the format precision for the channel_E0 value
        int_part_length = len(str(int(channel_E0)).replace(".", ""))
        format_precision = max(err_significant_digits - int_part_length, 0)
        format_str = f"{{:.{format_precision}f}}"
        channel_E0_str = format_str.format(channel_E0)

        # print(len(str(err_str)))
        # print(err_str)

        # Combine the channel_E0 value with the error part
        channel_E0_with_error = f"{channel_E0_str}({err_str})"
    else:
        channel_E0_with_error = "-"
    return channel_E0_with_error


# Read CSV files
metadata = pd.read_csv("./metadata/metadata_spectralDensity.csv")
# f_meson_gevp = pd.read_csv('./CSVs/F_meson_GEVP.csv')
# as_meson_gevp = pd.read_csv('./CSVs/AS_meson_GEVP.csv')

ensembles = ["M1", "M2", "M3", "M4", "M5"]
prefix = [
    "48x20x20x20b6.5mf0.71mas1.01",
    "64x20x20x20b6.5mf0.71mas1.01",
    "96x20x20x20b6.5mf0.71mas1.01",
    "64x20x20x20b6.5mf0.70mas1.01",
    "64x32x32x32b6.5mf0.72mas1.01",
]
prefix2 = [
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20",
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20",
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T96L20",
    "Sp4b6.5nF2nAS3mF-0.7mAS-1.01T64L20",
    "Sp4b6.5nF2nAS3mF-0.72mAS-1.01T64L32",
]

# Iterate through chunks of 4 rows in M3_spectral_density_spectrum.csv
chunk_size = 4


for n in range(3):
    if n == 0:
        energy = "E0"
    if n == 1:
        energy = "E1"
    if n == 2:
        energy = "E2"

    for index, ensemble in enumerate(ensembles):
        # Initialize LaTeX table string
        # latex_table = "\\begin{table}[ht]\n"
        # latex_table += "\\centering\n"
        latex_table = "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}\n"
        latex_table += "\\hline \\hline\n"
        latex_table += (
            "$C$ & $k$ & $N_{\\text{source}}$ & $N_{\\text{sink}}$ & "
            f"$aE_{n}$ $k$-G & $aE_{n}$ $(k+1)$-G & $aE_{n}$ $k$-C & $aE_{n}$ $(k+1)$-C$ & am_C & "
            "$\sigma_{G} / m_C$ & $\sigma_{C} / m_C$ \\\\\n"
        )
        latex_table += "\\hline\n"
        for chunk in pd.read_csv(
            f"./CSVs/{ensemble}_spectral_density_spectrum.csv", chunksize=chunk_size
        ):
            channel = chunk["channel"].min()
            # print(channel)
            repr = chunk["rep"].min()
            if repr == "fund":
                rep = "f"
            elif repr == "as":
                rep = repr

            ENSEMBLE = prefix2[index]
            CHANNEL2 = "ciao"
            if channel == "g5" and repr == "fund":
                CHANNEL2 = "ps"
                CHANNEL3 = "PS"
                try:
                    file_path = (
                        f"./JSONs/{ENSEMBLE}/meson_gevp_{rep}_{CHANNEL2}_samples.json"
                    )
                    print(file_path)
                    channel_E0, err_channel_E0 = process_gevp_En_mass_samples(
                        file_path, CHANNEL2, rep, n
                    )
                    # Check if the results are nan
                    if np.isnan(channel_E0) or np.isnan(err_channel_E0):
                        raise ValueError("Result contains NaN values")
                    err_channel_E0 = 30.0 * err_channel_E0
                    print(f"average: {channel_E0}")
                    print(f"error: {err_channel_E0}")
                except KeyError:
                    channel_E0 = "-"
                    err_channel_E0 = "-"
                except ValueError as e:
                    # Handle cases where the results are NaN
                    print(f"ValueError: {e}")
                    channel_E0 = "NaN"
                    err_channel_E0 = "NaN"

            # Add similar try-except blocks for other conditions
            elif channel == "g5" and repr == "as":
                CHANNEL2 = "ps"
                CHANNEL3 = CHANNEL2
                try:
                    file_path = (
                        f"./JSONs/{ENSEMBLE}/meson_gevp_{rep}_{CHANNEL2}_samples.json"
                    )
                    print(file_path)
                    channel_E0, err_channel_E0 = process_gevp_En_mass_samples(
                        file_path, CHANNEL2, rep, n
                    )
                    # Check if the results are nan
                    if np.isnan(channel_E0) or np.isnan(err_channel_E0):
                        raise ValueError("Result contains NaN values")
                    err_channel_E0 = 30.0 * err_channel_E0
                    print(f"average: {channel_E0}")
                    print(f"error: {err_channel_E0}")
                except KeyError:
                    channel_E0 = "-"
                    err_channel_E0 = "-"
                except ValueError as e:
                    # Handle cases where the results are NaN
                    print(f"ValueError: {e}")
                    channel_E0 = "NaN"
                    err_channel_E0 = "NaN"
                # Add similar try-except blocks for other conditions
            elif channel == "gi" and repr == "fund":
                CHANNEL2 = "v"
                CHANNEL3 = "V"
                try:
                    file_path = (
                        f"./JSONs/{ENSEMBLE}/meson_gevp_{rep}_{CHANNEL2}_samples.json"
                    )
                    print(file_path)
                    channel_E0, err_channel_E0 = process_gevp_En_mass_samples(
                        file_path, CHANNEL2, rep, n
                    )
                    # Check if the results are nan
                    if np.isnan(channel_E0) or np.isnan(err_channel_E0):
                        raise ValueError("Result contains NaN values")
                    err_channel_E0 = 30.0 * err_channel_E0
                    print(f"average: {channel_E0}")
                    print(f"error: {err_channel_E0}")
                except KeyError:
                    channel_E0 = "-"
                    err_channel_E0 = "-"
                except ValueError as e:
                    # Handle cases where the results are NaN
                    print(f"ValueError: {e}")
                    channel_E0 = "NaN"
                    err_channel_E0 = "NaN"
                # Add similar try-except blocks for other conditions
            elif channel == "gi" and repr == "as":
                CHANNEL2 = "v"
                CHANNEL3 = CHANNEL2
                try:
                    file_path = (
                        f"./JSONs/{ENSEMBLE}/meson_gevp_{rep}_{CHANNEL2}_samples.json"
                    )
                    print(file_path)
                    channel_E0, err_channel_E0 = process_gevp_En_mass_samples(
                        file_path, CHANNEL2, rep, n
                    )
                    # Check if the results are nan
                    if np.isnan(channel_E0) or np.isnan(err_channel_E0):
                        raise ValueError("Result contains NaN values")
                    err_channel_E0 = 30.0 * err_channel_E0
                    print(f"average: {channel_E0}")
                    print(f"error: {err_channel_E0}")
                except KeyError:
                    channel_E0 = "-"
                    err_channel_E0 = "-"
                except ValueError as e:
                    # Handle cases where the results are NaN
                    print(f"ValueError: {e}")
                    channel_E0 = "NaN"
                    err_channel_E0 = "NaN"
                # Add similar try-except blocks for other conditions
            elif channel == "g0gi" and repr == "fund":
                CHANNEL2 = "t"
                CHANNEL3 = "T"
                try:
                    file_path = (
                        f"./JSONs/{ENSEMBLE}/meson_gevp_{rep}_{CHANNEL2}_samples.json"
                    )
                    print(file_path)
                    channel_E0, err_channel_E0 = process_gevp_En_mass_samples(
                        file_path, CHANNEL2, rep, n
                    )
                    # Check if the results are nan
                    if np.isnan(channel_E0) or np.isnan(err_channel_E0):
                        raise ValueError("Result contains NaN values")
                    err_channel_E0 = 30.0 * err_channel_E0
                    print(f"average: {channel_E0}")
                    print(f"error: {err_channel_E0}")
                except KeyError:
                    channel_E0 = "-"
                    err_channel_E0 = "-"
                except ValueError as e:
                    # Handle cases where the results are NaN
                    print(f"ValueError: {e}")
                    channel_E0 = "NaN"
                    err_channel_E0 = "NaN"
                # Add similar try-except blocks for other conditions
            elif channel == "g0gi" and repr == "as":
                CHANNEL2 = "t"
                CHANNEL3 = CHANNEL2
                try:
                    file_path = (
                        f"./JSONs/{ENSEMBLE}/meson_gevp_{rep}_{CHANNEL2}_samples.json"
                    )
                    print(file_path)
                    channel_E0, err_channel_E0 = process_gevp_En_mass_samples(
                        file_path, CHANNEL2, rep, n
                    )
                    # Check if the results are nan
                    if np.isnan(channel_E0) or np.isnan(err_channel_E0):
                        raise ValueError("Result contains NaN values")
                    err_channel_E0 = 30.0 * err_channel_E0
                    print(f"average: {channel_E0}")
                    print(f"error: {err_channel_E0}")
                except KeyError:
                    channel_E0 = "-"
                    err_channel_E0 = "-"
                except ValueError as e:
                    # Handle cases where the results are NaN
                    print(f"ValueError: {e}")
                    channel_E0 = "NaN"
                    err_channel_E0 = "NaN"
            elif channel == "g5gi" and repr == "fund":
                CHANNEL2 = "av"
                CHANNEL3 = "AV"
                try:
                    file_path = (
                        f"./JSONs/{ENSEMBLE}/meson_gevp_{rep}_{CHANNEL2}_samples.json"
                    )
                    print(file_path)
                    channel_E0, err_channel_E0 = process_gevp_En_mass_samples(
                        file_path, CHANNEL2, rep, n
                    )
                    # Check if the results are nan
                    if np.isnan(channel_E0) or np.isnan(err_channel_E0):
                        raise ValueError("Result contains NaN values")
                    err_channel_E0 = 30.0 * err_channel_E0
                    print(f"average: {channel_E0}")
                    print(f"error: {err_channel_E0}")
                except KeyError:
                    channel_E0 = "-"
                    err_channel_E0 = "-"
                except ValueError as e:
                    # Handle cases where the results are NaN
                    print(f"ValueError: {e}")
                    channel_E0 = "NaN"
                    err_channel_E0 = "NaN"
                # Add similar try-except blocks for other conditions
            elif channel == "g5gi" and repr == "as":
                CHANNEL2 = "av"
                CHANNEL3 = CHANNEL2
                try:
                    file_path = (
                        f"./JSONs/{ENSEMBLE}/meson_gevp_{rep}_{CHANNEL2}_samples.json"
                    )
                    print(file_path)
                    channel_E0, err_channel_E0 = process_gevp_En_mass_samples(
                        file_path, CHANNEL2, rep, n
                    )
                    # Check if the results are nan
                    if np.isnan(channel_E0) or np.isnan(err_channel_E0):
                        raise ValueError("Result contains NaN values")
                    err_channel_E0 = 30.0 * err_channel_E0
                    print(f"average: {channel_E0}")
                    print(f"error: {err_channel_E0}")
                except KeyError:
                    channel_E0 = "-"
                    err_channel_E0 = "-"
                except ValueError as e:
                    # Handle cases where the results are NaN
                    print(f"ValueError: {e}")
                    channel_E0 = "NaN"
                    err_channel_E0 = "NaN"
                # Add similar try-except blocks for other conditions
            elif channel == "g0g5gi" and repr == "fund":
                CHANNEL2 = "at"
                CHANNEL3 = "AT"
                try:
                    file_path = (
                        f"./JSONs/{ENSEMBLE}/meson_gevp_{rep}_{CHANNEL2}_samples.json"
                    )
                    print(file_path)
                    channel_E0, err_channel_E0 = process_gevp_En_mass_samples(
                        file_path, CHANNEL2, rep, n
                    )
                    # Check if the results are nan
                    if np.isnan(channel_E0) or np.isnan(err_channel_E0):
                        raise ValueError("Result contains NaN values")
                    print(f"average: {channel_E0}")
                    err_channel_E0 = 30.0 * err_channel_E0
                    print(f"error: {err_channel_E0}")
                except KeyError:
                    channel_E0 = "-"
                    err_channel_E0 = "-"
                except ValueError as e:
                    # Handle cases where the results are NaN
                    print(f"ValueError: {e}")
                    channel_E0 = "NaN"
                    err_channel_E0 = "NaN"
                # Add similar try-except blocks for other conditions
            elif channel == "g0g5gi" and repr == "as":
                CHANNEL2 = "at"
                CHANNEL3 = CHANNEL2
                try:
                    file_path = (
                        f"./JSONs/{ENSEMBLE}/meson_gevp_{rep}_{CHANNEL2}_samples.json"
                    )
                    print(file_path)
                    channel_E0, err_channel_E0 = process_gevp_En_mass_samples(
                        file_path, CHANNEL2, rep, n
                    )
                    # Check if the results are nan
                    if np.isnan(channel_E0) or np.isnan(err_channel_E0):
                        raise ValueError("Result contains NaN values")
                    err_channel_E0 = 30.0 * err_channel_E0
                    print(f"average: {channel_E0}")
                    print(f"error: {err_channel_E0}")
                except KeyError:
                    channel_E0 = "-"
                    err_channel_E0 = "-"
                except ValueError as e:
                    # Handle cases where the results are NaN
                    print(f"ValueError: {e}")
                    channel_E0 = "NaN"
                    err_channel_E0 = "NaN"
                # Add similar try-except blocks for other conditions
            elif channel == "id" and repr == "fund":
                CHANNEL2 = "s"
                CHANNEL3 = "S"
                try:
                    file_path = (
                        f"./JSONs/{ENSEMBLE}/meson_gevp_{rep}_{CHANNEL2}_samples.json"
                    )
                    print(file_path)
                    channel_E0, err_channel_E0 = process_gevp_En_mass_samples(
                        file_path, CHANNEL2, rep, n
                    )
                    # Check if the results are nan
                    if np.isnan(channel_E0) or np.isnan(err_channel_E0):
                        raise ValueError("Result contains NaN values")
                    err_channel_E0 = 30.0 * err_channel_E0
                    print(f"average: {channel_E0}")
                    print(f"error: {err_channel_E0}")
                except KeyError:
                    channel_E0 = "-"
                    err_channel_E0 = "-"
                except ValueError as e:
                    # Handle cases where the results are NaN
                    print(f"ValueError: {e}")
                    channel_E0 = "NaN"
                    err_channel_E0 = "NaN"
                # Add similar try-except blocks for other conditions
            elif channel == "id" and repr == "as":
                CHANNEL2 = "s"
                CHANNEL3 = CHANNEL2
                try:
                    file_path = (
                        f"./JSONs/{ENSEMBLE}/meson_gevp_{rep}_{CHANNEL2}_samples.json"
                    )
                    print(file_path)
                    channel_E0, err_channel_E0 = process_gevp_En_mass_samples(
                        file_path, CHANNEL2, rep, n
                    )
                    # Check if the results are nan
                    if np.isnan(channel_E0) or np.isnan(err_channel_E0):
                        raise ValueError("Result contains NaN values")
                    err_channel_E0 = 30.0 * err_channel_E0
                    print(f"average: {channel_E0}")
                    print(f"error: {err_channel_E0}")
                except KeyError:
                    channel_E0 = "-"
                    err_channel_E0 = "-"
                except ValueError as e:
                    # Handle cases where the results are NaN
                    print(f"ValueError: {e}")
                    channel_E0 = "NaN"
                    err_channel_E0 = "NaN"
                # Add similar try-except blocks for other conditions

            # print(CHANNEL)
            # Extract required values from metadata
            k_peaks = metadata.loc[
                metadata["Ensemble"] == ensemble, f"{CHANNEL2}_k_peaks"
            ].values[0]
            n_source = metadata.loc[
                metadata["Ensemble"] == ensemble, f"{CHANNEL2}_Nsource_1"
            ].values[0]
            n_sink = metadata.loc[
                metadata["Ensemble"] == ensemble, f"{CHANNEL2}_Nsink_1"
            ].values[0]
            sigma1_over_m = metadata.loc[
                metadata["Ensemble"] == ensemble, f"{CHANNEL2}_sigma1_over_m"
            ].values[0]
            sigma2_over_m = metadata.loc[
                metadata["Ensemble"] == ensemble, f"{CHANNEL2}_sigma2_over_m"
            ].values[0]

            try:
                peak_gauss_min = chunk.loc[chunk["kernel"] == "GAUSS", "peaks"].min()
                peak_gauss_max = chunk.loc[chunk["kernel"] == "GAUSS", "peaks"].max()
                peak_cauchy_min = chunk.loc[chunk["kernel"] == "CAUCHY", "peaks"].min()
                peak_cauchy_max = chunk.loc[chunk["kernel"] == "CAUCHY", "peaks"].max()

                gauss_min = chunk.loc[
                    (chunk["kernel"] == "GAUSS") & (chunk["peaks"] == peak_gauss_min),
                    f"aE_{n}",
                ].min()
                err_gauss_min = chunk.loc[
                    (chunk["kernel"] == "GAUSS") & (chunk["peaks"] == peak_gauss_min),
                    f"errorE{n}",
                ].min()
                # print(gauss_min)
                gauss_max = chunk.loc[
                    (chunk["kernel"] == "GAUSS") & (chunk["peaks"] == peak_gauss_max),
                    f"aE_{n}",
                ].min()
                err_gauss_max = chunk.loc[
                    (chunk["kernel"] == "GAUSS") & (chunk["peaks"] == peak_gauss_max),
                    f"errorE{n}",
                ].min()
                # print(gauss_max)
                cauchy_min = chunk.loc[
                    (chunk["kernel"] == "CAUCHY") & (chunk["peaks"] == peak_cauchy_min),
                    f"aE_{n}",
                ].min()
                err_cauchy_min = chunk.loc[
                    (chunk["kernel"] == "CAUCHY") & (chunk["peaks"] == peak_cauchy_min),
                    f"errorE{n}",
                ].min()
                # print(cauchy_min)
                cauchy_max = chunk.loc[
                    (chunk["kernel"] == "CAUCHY") & (chunk["peaks"] == peak_cauchy_max),
                    f"aE_{n}",
                ].min()
                err_cauchy_max = chunk.loc[
                    (chunk["kernel"] == "CAUCHY") & (chunk["peaks"] == peak_cauchy_max),
                    f"errorE{n}",
                ].min()
                # print(cauchy_max, '\n\n')
                # print(channel_E0)
                if gauss_min != gauss_min:
                    gauss_min_with_error = "-"
                    gauss_max_with_error = "-"
                    cauchy_min_with_error = "-"
                    cauchy_max_with_error = "-"
                else:
                    if err_gauss_min != 0.0:
                        gauss_min_with_error = add_error(gauss_min, err_gauss_min)
                    else:
                        gauss_min_with_error = add_error(gauss_min, 0.01 * gauss_min)
                    if err_gauss_max != 0.0:
                        gauss_max_with_error = add_error(gauss_max, 0.01 * gauss_max)
                    if err_cauchy_min != 0.0:
                        cauchy_min_with_error = add_error(cauchy_min, err_cauchy_min)
                    else:
                        cauchy_min_with_error = add_error(cauchy_min, 0.01 * cauchy_min)
                    if err_cauchy_max != 0.0:
                        cauchy_max_with_error = add_error(cauchy_max, err_cauchy_max)
                    else:
                        cauchy_max_with_error = add_error(cauchy_max, 0.01 * cauchy_max)

                channel_E0_with_error = add_error(channel_E0, err_channel_E0)
            except KeyError:
                gauss_min_with_error = "-"
                gauss_max_with_error = "-"
                cauchy_min_with_error = "-"
                cauchy_max_with_error = "-"
                channel_E0_with_error = "-"

            # Adding the formatted value to the LaTeX table
            latex_table += f"{CHANNEL3} & {k_peaks} & {n_source} & {n_sink} & {gauss_min_with_error} & {gauss_max_with_error} & {cauchy_min_with_error} & {cauchy_max_with_error} & {channel_E0_with_error} & {sigma1_over_m} & {sigma2_over_m} \\\\\n"
            if CHANNEL2 == "s":
                latex_table += "\\hline \n"

        # Close LaTeX table for each chunk
        latex_table += "\\hline \\hline\n"
        latex_table += "\\end{tabular}\n"
        # latex_table += "\\caption{Your caption here.}\n"
        # latex_table += "\\label{table:my_table}\n"
        # latex_table += "\\end{table}\n"

        # Write LaTeX table to a file for each chunk
        with open(f"./assets/tables/{ensemble}_aE{n}_meson.tex", "w") as file:
            file.write(latex_table)
        # Reset LaTeX table for next chunk
        latex_table = ""

        # Print confirmation message
        print(f"Tables generated and saved in {ensemble}_aE{n}_meson.tex")
