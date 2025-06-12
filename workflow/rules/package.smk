from glob import glob

parsing_base = "plateaus/HiRep_parsing"
wall_parsing_base = "CB_autocorrelation_decay_constant/src_jl"
flow_parsing_base = "CB_autocorrelation_decay_constant/src_py"

topology_dirs = {
    "Lt48Ls20beta6.5mf0.71mas1.01FUN": "M1",
    "Lt64Ls20beta6.5mf0.71mas1.01FUN": "M2",
    "Lt96Ls20beta6.5mf0.71mas1.01FUN": "M3",
    "Lt64Ls20beta6.5mf0.70mas1.01FUN": "M4",
    "Lt64Ls32beta6.5mf0.72mas1.01FUN": "M5",
}

rule package_wall:
    input:
        metadata="metadata/ensemble_metadata.csv",
#         logs=...,
        instantiate_julia="intermediary_data/cb_julia_instantiated",
        script=f"{wall_parsing_base}/parse.jl",
    output:
        h5=protected("data_assets/wall_correlators.h5"),
    conda: "../envs/CB_autocorrelation_decay_constant.yml"
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


rule package_single_gflow:
    params:
        ensemble=lambda wildcards: topology_dirs[wildcards.subdir],
    input:
        log="raw_data/topology/{subdir}/out/out_flow",
        script=f"{flow_parsing_base}/package_flows_multirep.py",
    output:
        h5="intermediary_data/{subdir}.h5",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python {input.script} {input.log} --ensemble {params.ensemble} --h5_filename {output.h5}"


rule concatenate_gflow_hdf5:
    input:
        files=expand(
            "intermediary_data/{subdir}.h5",
            subdir=topology_dirs.keys(),
        ),
        script="src/hdf5_concatenate.py",
    output:
        h5=protected("data_assets/flows.h5"),
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python {input.script} {input.files} --output_file {output.h5}"
