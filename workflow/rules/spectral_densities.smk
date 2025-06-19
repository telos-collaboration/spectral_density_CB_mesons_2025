ensemble_prefixes = [
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20",
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20",
    "Sp4b6.5nF2nAS3mF-0.71mAS-1.01T96L20",
    "Sp4b6.5nF2nAS3mF-0.7mAS-1.01T64L20",
    "Sp4b6.5nF2nAS3mF-0.72mAS-1.01T64L32",
]
# chimera_channels = [
#     "Chimera_OC_even",
#     "Chimera_OC_odd",
#     "Chimera_OV12_even",
#     "Chimera_OV12_odd",
#     "Chimera_OV32_even",
#     "Chimera_OV32_odd",
# ]

chimera_channels = [
    "lambda_even",
    "sigma_even",
    "sigmastar_even",
    "lambda_odd",
    "sigma_odd",
    "sigmastar_odd",
]
meson_channels = ["ps", "v", "av", "t", "at", "s"]


rule analysis_template:
    threads: workflow.cores
    conda: "../envs/spectral_densities.yml"
    shell: "cd lsd_out && python ../{input.script}"


rule analysis_template_with_plot:
    threads: workflow.cores
    conda: "../envs/spectral_densities.yml"
    shell: "cd lsd_out && python ../{input.script} --plot_styles ../{input.plot_styles}"


use rule analysis_template as analyse_data_mesons with:
    input:
        script="lsd_out/analyse_data_mesons.py",
        data="input_correlators/chimera_data_reduced.h5",
        metadata="metadata/metadata_spectralDensity.csv",
    output:
        completion_tag="lsd_out/analyse_data_mesons_complete",
        useful_logs=expand(
            "input_fit/stability_plot/InverseProblemLOG_Alpha{index}.log",
            index=["A", "B", "C"],
        ),
        # samples1 = expand(
        #     "input_fit/lsdensity_samples1/lsdensitiesamplesE{energy:.16f}sig{sig:.16f}",
        #     zip,
        #     energy=[0.3, 0.465, 0.63, 0.795, 0.96, 1.125, 1.29],
        #     sig=[5.56e-7, 5.56e-7, 9.33e-7, 1.06e-6, 1.69e-6, 1.56e-6, 8.7e-7],
        # ),
        # samples2 = expand(
        #     "input_fit/lsdensity_samples2/lsdensitiesamplesE{energy:.16f}sig{sig:.16f}",
        #     zip,
        #     energy=[0.3, 0.465, 0.63, 0.795, 0.96, 1.125, 1.29],
        #     sig=[1.28e-6, 1.28e-6, 1.28e-6, 2.6e-6, 2.52e-6, 2.4e-6, 1.467e-6],
        # ),
    log: "lsd_out/paths.log"


use rule analysis_template as analyse_data_CB with:
    input:
        script="lsd_out/analyse_data_CB.py",
        data="input_correlators/chimera_data_reduced.h5",
        metadata="metadata/metadata_spectralDensity_chimerabaryons.csv",
    output:
        completion_tag="lsd_out/analyse_data_CB_complete",
    log: "lsd_out/paths.log"


use rule analysis_template as print_samples_mesons with:
    input:
        script="lsd_out/print_samples_mesons.py",
        metadata="metadata/metadata_spectralDensity.csv",
        dependency_tag="lsd_out/analyse_data_mesons_complete",
    output:
        completion_tag="lsd_out/print_samples_mesons_complete",
    log: "lsd_out/paths.log"


use rule analysis_template as print_samples_CB with:
    input:
        script="lsd_out/print_samples_CB.py",
        metadata="metadata/metadata_spectralDensity_chimerabaryons.csv",
        dependency_tag="lsd_out/analyse_data_CB_complete",
    output:
        completion_tag="lsd_out/print_samples_CB_complete",
    log: "lsd_out/paths.log"


rule fit_data_mesons:
    threads: 20
    input:
        script="lsd_out/fit_data_mesons.py",
        metadata_spectral_density="metadata/metadata_spectralDensity.csv",
        metadata_ratio_guesses="metadata/ratioguesses_spectrum.csv",
        dependency_tag="lsd_out/print_samples_mesons_complete",
    output:
        spectrum="CSVs/{ensemble,M[0-9]}_spectral_density_spectrum.csv",
    conda: "../envs/spectral_densities.yml"
    shell: "cd lsd_out && python ../{input.script} --ensembles {wildcards.ensemble}"


rule fit_data_CB:
    threads: 20
    input:
        script="lsd_out/fit_data_CB.py",
        metadata_spectral_density="metadata/metadata_spectralDensity.csv",
        metadata_ratio_guesses="metadata/ratioguesses_spectrum.csv",
        dependency_tag="lsd_out/print_samples_CB_complete",
    output:
        spectrum="CSVs/{ensemble}_chimerabaryons_spectral_density_spectrum.csv",
    conda: "../envs/spectral_densities.yml"
    shell: "cd lsd_out && python ../{input.script} --ensembles {wildcards.ensemble}"


use rule analysis_template_with_plot as simultaneous_fits_mesons with:
    input:
        script="lsd_out/simultaneous_fits_mesons.py",
        metadata="metadata/metadata_spectralDensity.csv",
        dependency_tag="lsd_out/fit_data_mesons_complete",
        plot_styles=plot_styles,
    output:
        completion_tag="lsd_out/simultaneous_fit_mesons_complete",
        matrix_elements=expand(
            "CSVs/{ensemble}_spectral_density_matrix_elements_CB.csv",
            ensemble=ensembles,
        ),


use rule analysis_template_with_plot as simultaneous_fits_CB with:
    input:
        script="lsd_out/simultaneous_fits_CB.py",
        metadata="metadata/metadata_spectralDensity_chimerabaryons.csv",
        dependency_tag="lsd_out/fit_data_CB_complete",
        plot_styles=plot_styles,
    output:
        completion_tag="lsd_out/simultaneous_fit_CB_complete",


rule post_analysis_spdens:
    input:
        script="lsd_out/post_analysis_spdens.py",
        metadata="metadata/renormalise.csv",
        mesons="lsd_out/simultaneous_fits_mesons_complete",
        cb="lsd_out/simultaneous_fits_CB_complete",
        topology="data_assets/flows.h5",
        matrix_elements=expand(
            "CSVs/{ensemble}_spectral_density_matrix_elements{group}.csv",
            ensemble=ensembles,
            group=["", "_CB"],
        ),
    output:
        spectrum=expand(
            "input_fit/final_spectrum/{group}{ensemble}_{state}.txt",
            group=["", "CB_"],
            ensemble=ensembles,
            state=["ground", "first", "second"]
        ),
        cb_matrixelement=expand(
            "input_fit/final_matrixel/{group}{ensemble}_ground.txt",
            group=["", "CB_"],
            ensemble=ensembles,
        ),


rule output_template:
    conda: "../envs/spectral_densities.yml"
    shell: "python {input.script} --plot_styles {input.plot_styles}"


use rule output_template as CSVs_to_tables_meson_GEVP with:
    input:
        script="src/CSVs_to_tables_meson_GEVP.py",
        csv=expand(
            "CSVs/{ensemble}_spectral_density_spectrum.csv",
            ensemble=ensembles,
        ),
        json=expand(
            "JSONs/{ensemble_prefix}/meson_gevp_{representation}_{channel}_samples.json",
            ensemble_prefix=ensemble_prefixes,
            representation=["f", "as"],
            channel=meson_channels,
        ),
    output:
        latex=expand(
            "assets/tables/{ensemble}_aE{n}_meson.tex",
            ensemble=ensembles,
            n=[0, 1, 2],
        ),


use rule output_template as CSVs_to_tables_CB_GEVP with:
    input:
        script="src/CSVs_to_tables_CB_GEVP.py",
        csvs=expand(
            "CSVs/{ensemble}_chimerabaryons_spectral_density_spectrum.csv",
            ensemble=ensembles,
        ),
        json=expand(
            "JSONs/{ensemble_prefix}/chimera_gevp_{channel}_samples.json",
            ensemble_prefix=ensemble_prefixes,
            channel=chimera_channels,
        ),
    output:
        latex=expand(
            "assets/tables/{ensemble}_aE{n}_CB.tex",
            ensemble=ensembles,
            n=[0, 1, 2],
        ),


use rule output_template as CSVs_to_tables_meson_matrixelements with:
    input:
        script="src/CSVs_to_tables_meson_matrixelements.py",
        matrixelement_csvs=expand(
            "CSVs/{ensemble}_spectral_density_matrix_elements.csv",
            ensemble=ensembles,
        ),
        spectrum_csvs=expand(
            "CSVs/{ensemble}_spectral_density_spectrum.csv",
            ensemble=ensembles,
        ),
    output:
        latex=expand(
            "assets/tables/{ensemble}_matrix_mesons.tex",
            ensemble=ensembles,
        )


use rule output_template as CSVs_to_tables_CB_matrixelements with:
    input:
        script="src/CSVs_to_tables_CB_matrixelements.py",
        metadata="metadata/metadata_spectralDensity_chimerabaryons.csv",
        csvs=expand(
            "CSVs/{ensemble}_spectral_density_matrix_elements_CB.csv",
            ensemble=ensembles,
        ),
        json=expand(
            "JSONs/{ensemble_prefix}/chimera_extraction_{channel}_samples.json",
            ensemble_prefix=ensemble_prefixes,
            channel=chimera_channels,
        )
    output:
        latex=expand(
            "assets/tables/{ensemble}_matrix_CB.tex",
            ensemble=ensembles,
        ),


rule renormalise:
    input:
        script="src/renormalise.py",
        data="metadata/renormalise.csv",
        mesons_tex=expand(
            "assets/tables/{ensemble}_matrix_mesons.tex",
            ensemble=ensembles,
        ),
        cb_tex=expand(
            "assets/tables/{ensemble}_matrix_CB.tex",
            ensemble=ensembles,
        ),
        topology="data_assets/flows.h5",
    output:
        mesons_tex=expand(
            "assets/tables/renormalised_{ensemble}_matrix_mesons.tex",
            ensemble=ensembles,
        ),
        cb_tex=expand(
            "assets/tables/renormalised_{ensemble}_matrix_CB.tex",
            ensemble=ensembles,
        ),
    conda: "../envs/spectral_densities.yml"
    shell: "python {input.script} --topology_h5 {input.topology}"


rule output_template_with_topology:
    conda: "../envs/spectral_densities.yml"
    shell: "python {input.script} --plot_styles {wildcards.plot_styles} --topology_h5 {input.topology}"


use rule output_template_with_topology as spectrum_MN_plot with:
    input:
        script="just_plotting/code/final_spectrum/spectrum_MN.py",
        ground_first_states=expand(
            "input_fit/final_spectrum/{state}{ensemble}_{excitation}.txt",
            state=["", "CB_"],
            ensemble=ensembles,
            excitation=["ground", "first"],
        ),
        second_states_meson=expand(
            "input_fit/final_spectrum/{ensemble}_second.txt",
            ensemble=["M2", "M3", "M4", "M5"],
        ),
        second_states_cb=expand(
            "input_fit/final_spectrum/CB_{ensemble}_second.txt",
            ensemble=ensembles,
        ),
        plot_styles=plot_styles,
        topology="data_assets/flows.h5",
    output:
        "assets/plots/final_spectrum_MN3.pdf",


use rule output_template_with_topology as spectrum_ensembles_showing_plot with:
    input:
        script="just_plotting/code/final_spectrum/spectrum_ensembles_showing.py",
        ground_first_states=expand(
            "input_fit/final_spectrum/{state}{ensemble}_{excitation}.txt",
            state=["", "CB_"],
            ensemble=ensembles,
            excitation=["ground", "first"],
        ),
        second_states_meson=expand(
            "input_fit/final_spectrum/{ensemble}_second.txt",
            ensemble=["M2", "M3", "M4", "M5"],
        ),
        second_states_cb=expand(
            "input_fit/final_spectrum/CB_{ensemble}_second.txt",
            ensemble=ensembles,
        ),
        plot_styles=plot_styles,
        topology="data_assets/flows.h5",
    output:
        "assets/plots/final_spectrum_detail.pdf",


use rule output_template_with_topology as matrix_MN_plot with:
    input:
        script="just_plotting/code/final_matrixel/matrix_MN.py",
        ground_states=expand(
            "input_fit/final_matrixel/{state}{ensemble}_ground.txt",
            state=["", "CB_"],
            ensemble=ensembles,
        ),
        plot_styles=plot_styles,
        topology="data_assets/flows.h5",
    output:
        "assets/plots/matrixel_fundamental_antisymmetric.pdf",
        "assets/plots/matrixel_chimera_baryons.pdf",


use rule output_template as stability_plot with:
    input:
        script="just_plotting/code/stability_plot/stability.py",
        data=expand(
            "input_fit/stability_plot/InverseProblemLOG_Alpha{index}.log",
            index=["A", "B", "C"],
        ),
        plot_styles=plot_styles,
    output:
        "assets/plots/LambdaScan.pdf",

use rule output_template as sfs_plot with:
    input:
        script="just_plotting/code/sfs/sfs.py",
        mesons="lsd_out/analyse_data_mesons_complete",
        cb="lsd_out/analyse_data_CB_complete",
        plot_styles=plot_styles,
    output:
        "assets/plots/spectral_density_single_gaussian.pdf",


rule weinberg:
    input:
        script="src/weinberg.py",
        metadata="metadata/renormalise.csv",
        json=expand(
            "JSONs/{ensemble_prefix}/meson_extraction_{representation}_{channel}_samples.json",
            ensemble_prefix=ensemble_prefixes,
            representation=["f", "as"],
            channel=["ps", "av", "v"],
        ),
        topology="data_assets/flows.h5",
    output:
        tex="assets/tables/s_parameters_table.tex",
        csv="CSVs/s_parameters_summary.csv",
    conda: "../envs/spectral_densities.yml"
    shell: "python {input.script} --topology_h5 {input.topology}"
