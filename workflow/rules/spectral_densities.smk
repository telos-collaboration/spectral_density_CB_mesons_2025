rule analysis_template:
    threads: workflow.cores
    conda: "../envs/spectral_densities.yml"
    shell: "cd lsd_out && python ../{input.script}"


use rule analysis_template as analyse_data_mesons with:
    input:
        script="lsd_out/analyse_data_mesons.py"
        data="input_correlators/chimera_data_reduced.hdf5",
        metadata="metadata/metadata_spectralDensity.csv",
    output:
        samples2 = expand(
            "input_fit/lsdensity_samples1/lsdensitiesamplesE{energy:.16f}sig{sig:.16f}",
            zip,
            energy=[0.3, 0.465, 0.63, 0.795, 0.96, 1.125, 1.29],
            sig=[5.56e-7, 5.56e-7, 9.33e-7, 1.06e-6, 1.69e-6, 1.56e-6, 8.7e-7],
        ),
        samples2 = expand(
            "input_fit/lsdensity_samples2/lsdensitiesamplesE{energy:.16f}sig{sig:.16f}",
            zip,
            energy=[0.3, 0.465, 0.63, 0.795, 0.96, 1.125, 1.29],
            sig=[1.28e-6, 1.28e-6, 1.28e-6, 2.6e-6, 2.52e-6, 2.4e-6, 1.467e-6],
        ),
    log: "lsd_out/paths.log"


use rule analysis_template as analyse_data_CB with:
    input:
        script="lsd_out/analyse_data_CB.py",
        data="input_correlators/chimera_data_reduced.hdf5",
        metadata="metadata/metadata_spectralDensity_chimerabaryons.csv",
    output: ...
    log: "lsd_out/paths.log"


rule output_template:
    conda: "../envs/spectral_densities.yml"
    shell: "python {input.script}"


use rule output_template as CSVs_to_tables_meson_GEVP with:
    input:
        script="CSVs_to_tables_meson_GEVP.py",
        data=...,
    output: ...


use rule output_template as CSVs_to_tables_CB_GEVP with:
    input:
        script="CSVs_to_tables_CB_GEVP.py",
        data=...,
    output: ...


use rule output_template as CSVs_to_tables_meson_GEVP with:
    input:
        script="CSVs_to_tables_meson_matrixelements.py",
        data=...,
    output: ...


use rule output_template as CSVs_to_tables_CB_GEVP with:
    input:
        script="CSVs_to_tables_CB_matrixelements.py",
        data=...,
    output: ...


use rule output_template as renormalise with:
    input:
        script="renormalise.py",
        data=...,
    output: ...


use rule output_template as spectrum_MN_plot with:
    input:
        script="just_plotting/code/final_spectrum/spectrum_MN.py",
        data=...,
    output: ...


use rule output_template as spectrum_ensembles_showing_plot with:
    input:
        script="just_plotting/code/final_spectrum/spectrum_ensembles_showing.py",
        data=...,
    output: ...


use rule output_template as matrix_MN_plot with:
    input:
        script="just_plotting/code/final_matrixel/matrix_MN.py",
        data=...,
    output: ...


use rule output_template as stability_plot with:
    input:
        script="just_plotting/code/final_matrixel/matrix_MN.py",
        data=...,
    output: ...


use rule output_template as elaborate_plot with:
    input:
        script="lsd_out/elaborate.py",
        data=...,
    output: ...
