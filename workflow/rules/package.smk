from glob import glob

parsing_base = "plateaus/HiRep_parsing"
wall_parsing_base = "CB_autocorrelation_decay_constant/src_jl"
flow_parsing_base = "CB_autocorrelation_decay_constant/src_py"


rule package_wall:
    input:
        metadata="metadata/ensembles.csv",
        logs=...,
        instantiate_julia="intermediary_data/cb_julia_instantiated",
        script=f"{wall_parsing_base}/parse.jl",
    output:
        h5: protected("data_assets/wall_correlators.h5"),
    conda: "../../CB_autocorrelation_decay_constant/environment.yml"
    shell: f"julia --project {wall_parsing_base} {{input.script}} --ensemble_metadata {{input.metadata}} --output_hdf5 {{output.h5}}"


rule package_smeared:
    params:
        file_dir="raw_data/corr/{smearing}",
        script_file_name="scripts/write_meson_{smearing}.jl"
    input:
        files=glob("raw_data/corr/{smearing}/*"),
        script=f"{parsing_base}/scripts/write_meson_{{smearing}}.jl",
    output:
        h5=protected("data_assets/correlators_{smearing}.h5"),
    conda:
        "../envs/hirep_parsing.yml"
    # Start packaging early,
    # since it is time consuming and many other processes depend on it
    priority:
        10
    shell:
        "cd {parsing_base} && julia {params.script_file_name} ../../{params.file_dir} ../../{output.h5}"


rule package_gflow:
    input:
        files=glob("raw_data/topology/*/out/out_flow"),
        script=f"{flow_parsing_base}/package_flows_multirep.py",
    output:
        h5=protected("data_assets/flows.h5"),
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python {input.script} {input.files} --h5_filename {output.h5}"
