import pandas as pd
import numpy as np

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
        channel_E0_with_error = '-'
    return channel_E0_with_error

# Read CSV files
metadata = pd.read_csv('./lsd_out/metadata/metadata_spectralDensity.csv')
F_meson_decay = pd.read_csv('./CSVs/F_meson_decay.csv')
AS_meson_decay = pd.read_csv('./CSVs/AS_meson_decay.csv')
ensembles = ['M1', 'M2', 'M3', 'M4', 'M5']

ensemble_map = {
    'M1': '48x20x20x20b6.5mf0.71mas1.01',
    'M2': '64x20x20x20b6.5mf0.71mas1.01',
    'M3': '96x20x20x20b6.5mf0.71mas1.01',
    'M4': '64x20x20x20b6.5mf0.70mas1.01',
    'M5': '64x32x32x32b6.5mf0.72mas1.01'
}

channel_map = {
    'g5': 'b_PS',
    'gi': 'b_V',
    'g0gi': 'b_T',
    'g5gi': 'b_AV',
    'g0g5gi': 'b_AT',
    'id': 'b_S'
}

channel_map2 = {
    'g5': 'PS',
    'gi': 'V',
    'g0gi': 'T',
    'g5gi': 'AV',
    'g0g5gi': 'AT',
    'id': 'S'
}

# Generate LaTeX table for each ensemble
for ensemble in ensembles:
    matrix_elements_path = f'./CSVs/{ensemble}_spectral_density_matrix_elements.csv'
    matrix_elements = pd.read_csv(matrix_elements_path)
    chunk_size = 4

    latex_table = "\\begin{table}[ht]\n"
    latex_table += "\\centering\n"
    latex_table += "\\begin{tabular}{|c|c|c|c|c|c|}\n"
    latex_table += "\\hline\n"
    latex_table += "$C$ & $ac_{n}$ G & $ac_{n}$ C & $ac_{0} (GEVP)$ & $\sigma_{G} / m_C$ & $\sigma_{C} / m_C$ \\\\\n"
    latex_table += "\\hline\n"

    row_count = 0

    for chunk in pd.read_csv(f'./CSVs/{ensemble}_spectral_density_spectrum.csv', chunksize=chunk_size):
        unique_channels = chunk['channel'].unique()

        for channel in unique_channels:
            CHANNEL2 = channel_map2.get(channel, 'Unknown')
            ch = CHANNEL2 if CHANNEL2 != 'Unknown' else channel
            if ch == 'PS' and row_count >= 6:
                ch = 'ps'
            if ch == 'V' and row_count >= 6:
                ch = 'v'
            if ch == 'T' and row_count >= 6:
                ch = 't'
            if ch == 'AV' and row_count >= 6:
                ch = 'av'
            if ch == 'AT' and row_count >= 6:
                ch = 'at'
            if ch == 'S' and row_count >= 6:
                ch = 's'
            #print('CHANNEL2', channel)
            # Safely retrieve metadata values
            try:
                sigma1_over_m = metadata.loc[metadata['Ensemble'] == ensemble, f"{ch}_sigma1_over_m"].values
                sigma2_over_m = metadata.loc[metadata['Ensemble'] == ensemble, f"{ch}_sigma2_over_m"].values
                sigma1_over_m = sigma1_over_m[0] if sigma1_over_m.size > 0 else '-'
                sigma2_over_m = sigma2_over_m[0] if sigma2_over_m.size > 0 else '-'
            except KeyError as e:
                print(f"KeyError for {ensemble} in metadata: {e}")
                sigma1_over_m = sigma2_over_m = '-'

            if row_count < 6:
                gauss_data = matrix_elements[(matrix_elements['kernel'] == 'GAUSS') & (matrix_elements['channel'] == channel)]['c0'].min()
                cauchy_data = matrix_elements[(matrix_elements['kernel'] == 'CAUCHY') & (matrix_elements['channel'] == channel)]['c0'].min()
            else:
                gauss_data = matrix_elements[(matrix_elements['kernel'] == 'GAUSS') & (matrix_elements['channel'] == channel)]['c0'].max()
                cauchy_data = matrix_elements[(matrix_elements['kernel'] == 'CAUCHY') & (matrix_elements['channel'] == channel)]['c0'].max()

            err_gauss_min = matrix_elements[(matrix_elements['kernel'] == 'GAUSS') & (matrix_elements['channel'] == channel)]['errorc0'].min()
            err_cauchy_min = matrix_elements[(matrix_elements['kernel'] == 'CAUCHY') & (matrix_elements['channel'] == channel)]['errorc0'].min()

            gauss_min_with_error = add_error(gauss_data, err_gauss_min) if not pd.isna(gauss_data) else '-'
            cauchy_min_with_error = add_error(cauchy_data, err_cauchy_min) if not pd.isna(cauchy_data) else '-'

            if row_count < 6:
                ac0_gevp_values = F_meson_decay.loc[
                    (F_meson_decay['ENS'] == ensemble_map[ensemble]),
                    [channel_map[channel], f"{channel_map[channel]}_error"]
                ].values
            else:
                ac0_gevp_values = AS_meson_decay.loc[
                    (AS_meson_decay['ENS'] == ensemble_map[ensemble]),
                    [channel_map[channel], f"{channel_map[channel]}_error"]
                ].values

            if ac0_gevp_values.size > 0:
                ac0_val, ac0_err = ac0_gevp_values[0]
                ac0_with_error = add_error(ac0_val, ac0_err)
            else:
                ac0_with_error = '-'

            latex_table += f"{ch} & {gauss_min_with_error} & {cauchy_min_with_error} & {ac0_with_error} & {sigma1_over_m} & {sigma2_over_m} \\\\\n"
            row_count += 1

    latex_table += "\\hline\n"
    latex_table += "\\end{tabular}\n"
    latex_table += "\\end{table}\n"

    with open(f'./tables/{ensemble}_output_table_matrix_mesons.tex', 'w') as file:
        file.write(latex_table)

    print(f"Table generated and saved in ./tables/{ensemble}_output_table_matrix_mesons.tex")
