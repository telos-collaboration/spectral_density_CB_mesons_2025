cb_julia_source = "CB_autocorrelation_decay_constant/src_jl"
cb_python_source = "CB_autocorrelation_decay_constant/src_py"

rule instantiate_julia:
    input:
        script=f"{cb_julia_source}/instantiate.jl",
    output:
        marker="intermediary_data/cb_julia_instantiated"
    conda: "../envs/CB_autocorrelation_decay_constant.yml"
    shell: f"cd {cb_julia_source}/.. && julia ../{{input.script}} && cd .. && touch {{output.marker}}"


rule analyse_flows:
    input:
        flows="data_assets/flows.h5",
        correlators="data_assets/wall_correlators.h5",
        instantiate_julia="intermediary_data/cb_julia_instantiated",
    output:
        table="assets/tables/ensembles.tex",
        h5="data_assets/autocorr.h5",
    conda: "../envs/CB_autocorrelation_decay_constant.yml"
    shell: f"julia --project {cb_julia_source} {cb_julia_source}/autocorrelation.jl --wilson_flow_hdf5 {{input.flows}} --wall_correlators_hdf5 {{input.correlators}} --output_hdf5 {{output.h5}} --output_tex {{output.table}} --plot_dir extra_assets/autocorrelation_plots"


rule mass_wall:
    params:
        plateau_start=lambda wildcards: lookup(within=metadata, query="ensemble_name == {wildcards.ensemble_name}", cols=f"wall_{wildcards.representation}_{wildcards.channel}_plateau_start"),
        plateau_end=lambda wildcards: lookup(within=metadata, query="ensemble_name == {wildcards.ensemble_name}", cols=f"wall_{wildcards.representation}_{wildcards.channel}_plateau_end"),
    input:
        correlators="data_assets/wall_correlators.h5",
        script=f"{cb_python_source}/mass_wall.py",
    output:
        csv="intermediary_data/wall_fits/{ensemble_name}{representation}_{channel}.csv",
    conda: "../envs/CB_autocorrelation_decay_constant.yml"
    shell: "python {input.script} --ensemble_name {wildcards.ensemble_name}/{wildcards.representation} --plateau_start {params.plateau_start} --plateau_end {params.plateau_end} --output_file_mean {output.csv} --channel {wildcards.channel} {input.correlators}"


rule comparison_plots:
    input:
        script=f"{cb_julia_source}/compare.jl",
        instantiate_julia="intermediary_data/cb_julia_instantiated",
        correlators="data_assets/wall_correlators.h5",
    output:
        comparison_csv="data_assets/comparison_table.csv",
        comparison_table="assets/tables/local_smeared_decay_constants.tex",
        comparison_plot="assets/plots/wall_comparison.pdf",
    conda: "../envs/CB_autocorrelation_decay_constant.yml"
    shell:
        f"julia --project {cb_julia_source} {{input.script}} --wall_correlators_h5 {{input.correlators}} --wall_fits intermediary_data/wall_fits --smeared_results external_data/smeared --decay_output_csv {{output.comparison_csv}} --decay_output_tex {{output.comparison_table}} --decay_output_pdf {{output.wall_comparison.pdf}}"
