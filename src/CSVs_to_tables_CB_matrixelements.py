import pandas as pd
import numpy as np
import json

def bootstrap_error(data, n_bootstrap=1000):
    """Compute the bootstrap error for a given dataset."""
    return np.std(data)


def add_error(channel_E0, err):
    if channel_E0 != 0 and not np.isnan(channel_E0):
        err_decimal_places = max(0, -int(np.floor(np.log10(err)))) if err > 0 else 0

        if err_decimal_places > 0:
            decimal_format = f".6f"
        else:
            decimal_format = ".0f"

        channel_E0_str = f"{channel_E0:{decimal_format}}"
        err_int = int(round(err * 1e6))

        channel_E0_with_error = f"{channel_E0_str}({err_int})"
    else:
        channel_E0_with_error = "-"
    return channel_E0_with_error


# Read metadata CSV
metadata = pd.read_csv("./metadata/metadata_spectralDensity_chimerabaryons.csv")
ensembles = ["M1", "M2", "M3", "M4", "M5"]
# ensembles = ['M1']
# Mapping for ensemble and channel names
ensemble_map = {
    "M1": "48x20x20x20b6.5mf0.71mas1.01",
    "M2": "64x20x20x20b6.5mf0.71mas1.01",
    "M3": "96x20x20x20b6.5mf0.71mas1.01",
    "M4": "64x20x20x20b6.5mf0.70mas1.01",
    "M5": "64x32x32x32b6.5mf0.72mas1.01",
}

channel_map = {
    "Chimera_OC_even": "b_Lambda_even",
    "Chimera_OV12_even": "b_Sigma_even",
    "Chimera_OV32_even": "b_SigmaS_even",
    "Chimera_OC_odd": "b_Lambda_odd",
    "Chimera_OV12_odd": "b_Sigma_odd",
    "Chimera_OV32_odd": "b_SigmaS_odd",
}

channel2 = [
    "lambda_even",
    "sigma_even",
    "sigmastar_even",
    "lambda_odd",
    "sigma_odd",
    "sigmastar_odd",
]
# channel2 = ['lambda_even', 'lambda_odd', 'sigma_even', 'sigma_odd', 'sigmastar_even', 'sigmastar_odd']

prefix = [
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20",
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20",
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T96L20",
    "Sp4b6.5nF2nAS3mF-0.7mAS-1.01T64L20",
    "Sp4b6.5nF2nAS3mF-0.72mAS-1.01T64L32",
]


# volumes = [20**3, 20**3, 20**3, 20**3, 32**3]


# Load the metadata
df = pd.read_csv("./metadata/ensemble_metadata.csv")

# Compute volume dictionary for ensembles M1–M5
volumes = {
    row["ensemble_name"]: int(row["Ns"]) ** 3
    for _, row in df.iterrows()
    if row["ensemble_name"] in ["M1", "M2", "M3", "M4", "M5"]
}


# Iterate through each ensemble
for idx, ensemble in enumerate(ensembles):
    ensemble2 = prefix[idx]

    matrix_elements_path = f"./CSVs/{ensemble}_spectral_density_matrix_elements_CB.csv"
    matrix_elements = pd.read_csv(matrix_elements_path)
    chunk_size = 2

    # Initialize LaTeX table string
    # latex_table = "\\begin{table}[ht]\n"
    # latex_table += "\\centering\n"
    latex_table = "\\begin{tabular}{|c|c|c|c|c|c|}\n"
    latex_table += "\\hline \\hline\n"
    latex_table += "$C$ & $a^{3}K_{B, 0}$-G & $a^{3}K_{B,0}$-C & $a^{3}K_{B,0} (\hbox{corr.})$ & $\sigma_{G} / m_C$ & $\sigma_{C} / m_C$ \\\\\n"
    latex_table += "\\hline\n"
    idx2 = 0
    # Read data for the current ensemble in chunks
    for chunk in pd.read_csv(
        f"./CSVs/{ensemble}_spectral_density_matrix_elements_CB.csv",
        chunksize=chunk_size,
    ):
        unique_channels = chunk["channel"].unique()

        for channel in unique_channels:
            # print(idx2)
            channelone = channel2[idx2]
            # Map channel to label and metadata keys
            if channel == "Chimera_OC_even":
                CHANNEL2, ch = "PS", "$\\Lambda^{+}_{\\rm CB}$"
            elif channel == "Chimera_OC_odd":
                CHANNEL2, ch = "V", "$\\Lambda^{-}_{\\rm CB}$"
            elif channel == "Chimera_OV12_even":
                CHANNEL2, ch = "T", "$\\Sigma^{+}_{\\rm CB}$"
            elif channel == "Chimera_OV12_odd":
                CHANNEL2, ch = "AV", "$\\Sigma^{-}_{\\rm CB}$"
            elif channel == "Chimera_OV32_even":
                CHANNEL2, ch = "AT", "$\\Sigma^{* \\, +}_{\\rm CB}$"
            elif channel == "Chimera_OV32_odd":
                CHANNEL2, ch = "S", "$\\Sigma^{* \\, -}_{\\rm CB}$"
            # Load JSON data
            # print(channelone)
            with open(
                f"./JSONs/{ensemble2}/chimera_extraction_{channelone}_samples.json", "r"
            ) as json_file:
                json_data = json.load(json_file)
            vol = volumes[ensemble]
            # Retrieve metadata values for the current channel and ensemble
            try:
                sigma1_over_m = metadata.loc[
                    metadata["Ensemble"] == ensemble, f"{CHANNEL2}_sigma1_over_m"
                ].values[0]
                sigma2_over_m = metadata.loc[
                    metadata["Ensemble"] == ensemble, f"{CHANNEL2}_sigma2_over_m"
                ].values[0]
            except KeyError as e:
                print(f"KeyError for {ensemble} in metadata: {e}")
                sigma1_over_m = sigma2_over_m = "-"

            # Retrieve c0 and errorc0 values for current channel in matrix elements
            gauss_data = matrix_elements[
                (matrix_elements["kernel"] == "GAUSS")
                & (matrix_elements["channel"] == channel)
            ]
            cauchy_data = matrix_elements[
                (matrix_elements["kernel"] == "CAUCHY")
                & (matrix_elements["channel"] == channel)
            ]

            # Check if data exists for GAUSS and CAUCHY kernels
            if not gauss_data.empty:
                gauss_min, err_gauss_min = (
                    gauss_data["c0"].min() / np.sqrt(vol),
                    gauss_data["errorc0"].min() / np.sqrt(vol),
                )
                gauss_min_with_error = add_error(gauss_min, err_gauss_min)
            else:
                gauss_min_with_error = "-"

            if not cauchy_data.empty:
                cauchy_min, err_cauchy_min = (
                    cauchy_data["c0"].min() / np.sqrt(vol),
                    cauchy_data["errorc0"].min() / np.sqrt(vol),
                )
                cauchy_min_with_error = add_error(cauchy_min, err_cauchy_min)
            else:
                cauchy_min_with_error = "-"
            print(f"{channelone}_matrix_element")
            if f"{channelone}_matrix_element" in json_data:
                samples = np.array(json_data[f"{channelone}_matrix_element"])
                ac0_val = np.mean(samples) / np.sqrt(vol)
                ac0_err = bootstrap_error(samples) / np.sqrt(vol)
                ac0_with_error = add_error(ac0_val, ac0_err)
            else:
                ac0_with_error = "-"
            idx2 = idx2 + 1
            # print('vol: ', vol)
            # Add the unique row for each channel to the LaTeX table
            latex_table += f"{ch} & {gauss_min_with_error} & {cauchy_min_with_error} & {ac0_with_error} & {sigma1_over_m} & {sigma2_over_m} \\\\\n"

    # Finalize LaTeX table and write to file
    latex_table += "\\hline \\hline\n"
    latex_table += "\\end{tabular}\n"
    # latex_table += "\\end{table}\n"

    with open(f"./assets/tables/{ensemble}_matrix_CB.tex", "w") as file:
        file.write(latex_table)

    print(f"Table generated and saved in {ensemble}_matrix_CB.tex")
