import pandas as pd
import numpy as np

def add_error(channel_E0, err):
    if channel_E0 != 0 and not np.isnan(channel_E0):
        # Calculate the number of decimal places in the error
        err_decimal_places = max(0, -int(np.floor(np.log10(err)))) if err > 0 else 0
        
        # Format the main value with the right number of decimal places
        if err_decimal_places > 0:
            decimal_format = f".{err_decimal_places}f"  # Show one more decimal place than the error
        else:
            decimal_format = ".0f"  # Show no decimals if error is very small

        channel_E0_str = f"{channel_E0:{decimal_format}}"  # Main value formatting

        # Round the error to the nearest integer for display
        err_int = int(round(err * 1e6))  # Convert to an integer in the same precision scale

        channel_E0_with_error = f"{channel_E0_str}({err_int})"

    else:
        channel_E0_with_error = '-'
    return channel_E0_with_error



# Read CSV files
metadata = pd.read_csv('./lsd_out/metadata/metadata_spectralDensity_chimerabaryons.csv')
CB_matrix_element = pd.read_csv('./CSVs/CB_matrix_element.csv')
ensembles = ['M1', 'M2', 'M3', 'M4', 'M5']

# Mapping for ensemble and channel names
ensemble_map = {
    'M1': '48x20x20x20b6.5mf0.71mas1.01',
    'M2': '64x20x20x20b6.5mf0.71mas1.01',
    'M3': '96x20x20x20b6.5mf0.71mas1.01',
    'M4': '64x20x20x20b6.5mf0.70mas1.01',
    'M5': '64x32x32x32b6.5mf0.72mas1.01'
}

channel_map = {
    'Chimera_OC_even': 'b_Lambda_even',
    'Chimera_OC_odd': 'b_Lambda_odd',
    'Chimera_OV12_even': 'b_Sigma_even',
    'Chimera_OV12_odd': 'b_Sigma_odd',
    'Chimera_OV32_even': 'b_SigmaS_even',
    'Chimera_OV32_odd': 'b_SigmaS_odd'
}

# Iterate through each ensemble
for ensemble in ensembles:
    matrix_elements_path = f'./CSVs/{ensemble}_spectral_density_matrix_elements_CB.csv'
    matrix_elements = pd.read_csv(matrix_elements_path)
    chunk_size = 4

    # Initialize LaTeX table string
    latex_table = "\\begin{table}[ht]\n"
    latex_table += "\\centering\n"
    latex_table += "\\begin{tabular}{|c|c|c|c|c|c|}\n"
    latex_table += "\\hline\n"
    latex_table += "$C$ & $ac_{n}$ G & $ac_{n}$ C & $ac_{0} (GEVP)$ & $\sigma_{G} / m_C$ & $\sigma_{C} / m_C$ \\\\\n"
    latex_table += "\\hline\n"

    # Read data for the current ensemble in chunks
    for chunk in pd.read_csv(f'./CSVs/{ensemble}_chimerabaryons_spectral_density_spectrum.csv', chunksize=chunk_size):
        unique_channels = chunk['channel'].unique()

        for channel in unique_channels:
            # Map channel to label and metadata keys
            if channel == 'Chimera_OC_even':
                CHANNEL2, ch = 'PS', '$\Lambda^{+}_{\\rm CB}$'
            elif channel == 'Chimera_OC_odd':
                CHANNEL2, ch = 'V', '$\Lambda^{-}_{\\rm CB}$'
            elif channel == 'Chimera_OV12_even':
                CHANNEL2, ch = 'T', '$\Sigma^{+}_{\\rm CB}$'
            elif channel == 'Chimera_OV12_odd':
                CHANNEL2, ch = 'AV', '$\Sigma^{-}_{\\rm CB}$'
            elif channel == 'Chimera_OV32_even':
                CHANNEL2, ch = 'AT', '$\Sigma^{* \, +}_{\\rm CB}$'
            elif channel == 'Chimera_OV32_odd':
                CHANNEL2, ch = 'S', '$\Sigma^{* \, -}_{\\rm CB}$'

            # Retrieve metadata values for the current channel and ensemble
            try:
                sigma1_over_m = metadata.loc[metadata['Ensemble'] == ensemble, f"{CHANNEL2}_sigma1_over_m"].values[0]
                sigma2_over_m = metadata.loc[metadata['Ensemble'] == ensemble, f"{CHANNEL2}_sigma2_over_m"].values[0]
            except KeyError as e:
                print(f"KeyError for {ensemble} in metadata: {e}")
                sigma1_over_m = sigma2_over_m = '-'

            # Retrieve c0 and errorc0 values for current channel in matrix elements
            gauss_data = matrix_elements[(matrix_elements['kernel'] == 'GAUSS') & (matrix_elements['channel'] == channel)]
            cauchy_data = matrix_elements[(matrix_elements['kernel'] == 'CAUCHY') & (matrix_elements['channel'] == channel)]
            
            # Check if data exists for GAUSS and CAUCHY kernels
            if not gauss_data.empty:
                gauss_min, err_gauss_min = gauss_data['c0'].min(), gauss_data['errorc0'].min()
                gauss_min_with_error = add_error(gauss_min, err_gauss_min)
            else:
                gauss_min_with_error = '-'

            if not cauchy_data.empty:
                cauchy_min, err_cauchy_min = cauchy_data['c0'].min(), cauchy_data['errorc0'].min()
                cauchy_min_with_error = add_error(cauchy_min, err_cauchy_min)
            else:
                cauchy_min_with_error = '-'

            # Retrieve the b_* values from CB_matrix_element.csv
            ac0_gevp_values = CB_matrix_element.loc[
                (CB_matrix_element['ENS'] == ensemble_map[ensemble]),
                [channel_map[channel], f"{channel_map[channel]}_error"]
            ].values

            # Format ac0 (GEVP) value with error
            if ac0_gevp_values.size > 0:
                ac0_val, ac0_err = ac0_gevp_values[0]
                ac0_with_error = add_error(ac0_val, ac0_err)
            else:
                ac0_with_error = '-'

            # Add the unique row for each channel to the LaTeX table
            latex_table += f"{ch} & {gauss_min_with_error} & {cauchy_min_with_error} & {ac0_with_error} & {sigma1_over_m} & {sigma2_over_m} \\\\\n"

    # Finalize LaTeX table and write to file
    latex_table += "\\hline\n"
    latex_table += "\\end{tabular}\n"
    latex_table += "\\end{table}\n"

    with open(f'./tables/{ensemble}_output_table_matrix_CB.tex', 'w') as file:
        file.write(latex_table)

    print(f"Table generated and saved in ./tables/{ensemble}_output_table_matrix_CB.tex")

