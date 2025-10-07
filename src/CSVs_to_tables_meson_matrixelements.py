import pandas as pd
import numpy as np
import json
import math

# TODO: Figure out if this is incorrect?
def bootstrap_error(data, n_resamples=1000):
    return np.std(data)


def add_error(channel_E0, err):
    if channel_E0 != 0 and not np.isnan(channel_E0):
        err_decimal_places = max(0, -int(np.floor(np.log10(err)))) if err > 0 else 0

        if err_decimal_places > 0:
            decimal_format = f".4f"
        else:
            decimal_format = ".0f"

        channel_E0_str = f"{channel_E0:{decimal_format}}"
        err_int = int(round(err * 1e4))

        channel_E0_with_error = f"{channel_E0_str}({err_int})"
    else:
        channel_E0_with_error = "-"
    return channel_E0_with_error


# Read metadata CSV file
metadata = pd.read_csv("./metadata/metadata_spectralDensity.csv")

ensemble_map = {
    "M1": "48x20x20x20b6.5mf0.71mas1.01",
    "M2": "64x20x20x20b6.5mf0.71mas1.01",
    "M3": "96x20x20x20b6.5mf0.71mas1.01",
    "M4": "64x20x20x20b6.5mf0.70mas1.01",
    "M5": "64x32x32x32b6.5mf0.72mas1.01",
}

channel_map = {
    "g5": "b_PS",
    "gi": "b_V",
    "g0gi": "b_T",
    "g5gi": "b_AV",
    "g0g5gi": "b_AT",
    "id": "b_S",
}

channel_map2 = {
    "g5": "PS",
    "gi": "V",
    "g0gi": "T",
    "g5gi": "AV",
    "g0g5gi": "AT",
    "id": "S",
}

# Generate LaTeX table for each ensemble
ensembles = ["M1", "M2", "M3", "M4", "M5"]
prefix = [
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20",
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20",
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T96L20",
    "Sp4b6.5nF2nAS3mF-0.7mAS-1.01T64L20",
    "Sp4b6.5nF2nAS3mF-0.72mAS-1.01T64L32",
]
for idx, ensemble in enumerate(ensembles):
    ensemble2 = prefix[idx]

    matrix_elements_path = f"./CSVs/{ensemble}_spectral_density_matrix_elements.csv"
    matrix_elements = pd.read_csv(matrix_elements_path)
    chunk_size = 4

    # latex_table = "\\begin{table}[ht]\n"
    # latex_table += "\\centering\n"
    latex_table = "\\begin{tabular}{|c|c|c|c|c|c|}\n"
    latex_table += "\\hline \\hline\n"
    latex_table += "$C$ & $a^{2}c_{M, 0}$-G & $a^{2}c_{M, 0}$-C & $a^{2}c_{M,0} (\hbox{corr.})$ & $\sigma_{G} / m_C$ & $\sigma_{C} / m_C$ \\\\\n"
    latex_table += "\\hline\n"

    row_count = 0
    rep = "fund"
    for chunk in pd.read_csv(
        f"./CSVs/{ensemble}_spectral_density_spectrum.csv", chunksize=chunk_size
    ):
        unique_channels = chunk["channel"].unique()
        print(unique_channels)
        for channel in unique_channels:
            CHANNEL2 = channel_map2.get(channel, "Unknown")
            ch = CHANNEL2 if CHANNEL2 != "Unknown" else channel
            if row_count == 6:
                row_count = 0
                rep = "as"
            if ch == "PS" and row_count >= 6:
                ch = "ps"
            if ch == "PS":
                ch2 = "ps"
            if ch == "V" and row_count >= 6:
                ch = "v"
            if ch == "V":
                ch2 = "v"
            if ch == "T" and row_count >= 6:
                ch = "t"
            if ch == "T":
                ch2 = "t"
            if ch == "AV" and row_count >= 6:
                ch = "av"
            if ch == "AV":
                ch2 = "av"
            if ch == "AT" and row_count >= 6:
                ch = "at"
            if ch == "AT":
                ch2 = "at"
            if ch == "S" and row_count >= 6:
                ch = "s"
            if ch == "S":
                ch2 = "s"

            # Safely retrieve metadata values
            try:
                sigma1_over_m = metadata.loc[
                    metadata["Ensemble"] == ensemble, f"{ch}_sigma1_over_m"
                ].values
                sigma2_over_m = metadata.loc[
                    metadata["Ensemble"] == ensemble, f"{ch}_sigma2_over_m"
                ].values
                sigma1_over_m = sigma1_over_m[0] if sigma1_over_m.size > 0 else "-"
                sigma2_over_m = sigma2_over_m[0] if sigma2_over_m.size > 0 else "-"
            except KeyError as e:
                print(f"KeyError for {ensemble} in metadata: {e}")
                sigma1_over_m = sigma2_over_m = "-"

                # Load JSON files
            with open(
                f"./JSONs/{ensemble2}/meson_extraction_f_{ch2}_samples.json", "r"
            ) as f:
                f_ps_data = json.load(f)

            with open(
                f"./JSONs/{ensemble2}/meson_extraction_as_{ch2}_samples.json", "r"
            ) as f:
                as_ps_data = json.load(f)

            if rep == "fund":
                print("ch2: ", ch2)
                name = f"f_{ch2}_matrix_element"
                data = f_ps_data[name]
            else:
                print("ch2: ", ch2)
                name = f"as_{ch2}_matrix_element"
                data = as_ps_data[name]

            if data:
                ac0_val = np.mean(data)
                ac0_err = bootstrap_error(data)
            else:
                ac0_val = ac0_err = 0

            gauss_data = matrix_elements[
                (matrix_elements["kernel"] == "GAUSS")
                & (matrix_elements["channel"] == channel)
                & (matrix_elements["rep"] == rep)
            ]["c0"].min()
            cauchy_data = matrix_elements[
                (matrix_elements["kernel"] == "CAUCHY")
                & (matrix_elements["channel"] == channel)
                & (matrix_elements["rep"] == rep)
            ]["c0"].min()

            err_gauss_min = matrix_elements[
                (matrix_elements["kernel"] == "GAUSS")
                & (matrix_elements["channel"] == channel)
                & (matrix_elements["rep"] == rep)
            ]["errorc0"].min()
            err_cauchy_min = matrix_elements[
                (matrix_elements["kernel"] == "CAUCHY")
                & (matrix_elements["channel"] == channel)
                & (matrix_elements["rep"] == rep)
            ]["errorc0"].min()
            print(rep)
            gauss_min_with_error = (
                add_error(gauss_data, err_gauss_min) if not pd.isna(gauss_data) else "-"
            )
            cauchy_min_with_error = (
                add_error(cauchy_data, err_cauchy_min)
                if not pd.isna(cauchy_data)
                else "-"
            )
            row_count += 1
            ac0_with_error = add_error(ac0_val, ac0_err)
            if rep == "fund":
                # Skip unwanted channels

                if ch in ["T", "AT", "S", "t", "at", "s"]:
                    continue
                latex_table += f"{ch} & {gauss_min_with_error} & {cauchy_min_with_error} & {ac0_with_error} & {sigma1_over_m} & {sigma2_over_m} \\\\ \n"
                if ch == "AV":
                    latex_table += "\\hline \n"
            else:
                # Skip unwanted channels
                if ch in ["T", "AT", "S", "t", "at", "s"]:
                    continue
                latex_table += f"{ch2} & {gauss_min_with_error} & {cauchy_min_with_error} & {ac0_with_error} & {sigma1_over_m} & {sigma2_over_m} \\\\ \n"

    latex_table += "\\hline \\hline\n"
    latex_table += "\\end{tabular}\n"
    # latex_table += "\\end{table}\n"

    with open(f"./assets/tables/{ensemble}_matrix_meson.tex", "w") as file:
        file.write(latex_table)

    print(f"Table generated and saved in ./assets/tables/{ensemble}_matrix_mesons.tex")
