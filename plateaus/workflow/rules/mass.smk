import pandas as pd
from functools import partial


all_metadata = pd.read_csv("../input_fit/metadata/ensemble_metadata.csv")

metadata_query = "Nc == {Nc} & Nt == {Nt} & Ns == {Ns} & beta == {beta} & nF == {nF} & mF == {mF} & nAS == {nAS} & mAS == {mAS}"
metadata_lookup = partial(lookup, within=all_metadata, query=metadata_query)

dir_template = "Sp{Nc}b{beta}nF{nF}nAS{nAS}mF{mF}mAS{mAS}T{Nt}L{Ns}"


rule gevp_chimera_baryon_mass:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=metadata_lookup(),
        E0_plateau_start=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E0_plateau_start"),
        E0_plateau_end=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E0_plateau_end"),
        E1_plateau_start=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E1_plateau_start"),
        E1_plateau_end=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E1_plateau_end"),
        E2_plateau_start=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E2_plateau_start"),
        E2_plateau_end=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E2_plateau_end"),
    input:
        data="../input_correlators/chimera_data_reduced.hdf5",
        script="src/mass_gevp_chimera.py",
    output:
        samples=f"../JSONs/{dir_template}/chimera_gevp_{{channel}}_samples.json",
        #plot=f"../JSONs/{dir_template}/chimera_gevp_{{channel}}_eff_mass.pdf",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_samples {output.samples} --ensemble_name {params.metadata.ensemble_name}"
        " --beta {params.metadata.beta} --mF {params.metadata.mF} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns}"
        " --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf_spectrum}"
        " --channel {wildcards.channel} --gevp_t0 {params.metadata.gevp_t0}"
        " --E0_plateau_start {params.E0_plateau_start} --E0_plateau_end {params.E0_plateau_end}"
        " --E1_plateau_start {params.E1_plateau_start} --E1_plateau_end {params.E1_plateau_end}"
        " --E2_plateau_start {params.E2_plateau_start} --E2_plateau_end {params.E2_plateau_end}"


rule chimera_matrix_element:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=metadata_lookup(),
        plateau_start=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_matrix_element_plateau_start"),
        plateau_end=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_matrix_element_plateau_end"),
    input:
        data="../input_correlators/chimera_data_reduced.hdf5",
        script="src/matrix_element_chimera.py",
    output:
        samples=f"../JSONs/{dir_template}/chimera_extraction_{{channel}}_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_samples {output.samples} --ensemble_name {params.metadata.ensemble_name}"
        " --beta {params.metadata.beta} --mF {params.metadata.mF} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns}"
        " --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf_spectrum}"
        " --channel {wildcards.channel} --E0_plateau_start {params.plateau_start} --E0_plateau_end {params.plateau_end}"


rule gevp_meson_mass:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=metadata_lookup(),
        E0_plateau_start=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E0_plateau_start"),
        E0_plateau_end=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E0_plateau_end"),
        E1_plateau_start=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E1_plateau_start"),
        E1_plateau_end=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E1_plateau_end"),
        E2_plateau_start=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E2_plateau_start"),
        E2_plateau_end=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_E2_plateau_end"),
    input:
        data="../input_correlators/chimera_data_reduced.hdf5",
        script="src/mass_gevp_meson.py",
    output:
        samples=f"../JSONs/{dir_template}/meson_gevp_{{channel}}_samples.json",
        #plot=f"../JSONs/{dir_template}/meson_gevp_{{channel}}_eff_mass.pdf",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_samples {output.samples} --ensemble_name {params.metadata.ensemble_name}"
        " --beta {params.metadata.beta} --mF {params.metadata.mF} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns}"
        " --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf_spectrum}"
        " --channel {wildcards.channel} --gevp_t0 {params.metadata.gevp_t0}"
        " --E0_plateau_start {params.E0_plateau_start} --E0_plateau_end {params.E0_plateau_end}"
        " --E1_plateau_start {params.E1_plateau_start} --E1_plateau_end {params.E1_plateau_end}"
        " --E2_plateau_start {params.E2_plateau_start} --E2_plateau_end {params.E2_plateau_end}"


rule meson_matrix_element:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
        metadata=metadata_lookup(),
        plateau_start=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_matrix_element_plateau_start"),
        plateau_end=lambda wildcards: metadata_lookup(cols=f"{wildcards.channel}_matrix_element_plateau_end"),
    input:
        data="../input_correlators/chimera_data_reduced.hdf5",
        script="src/matrix_element_meson.py",
    output:
        samples=f"../JSONs/{dir_template}/meson_extraction_{{channel}}_samples.json",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.data} --output_file_samples {output.samples} --ensemble_name {params.metadata.ensemble_name}"
        " --beta {params.metadata.beta} --mF {params.metadata.mF} --mAS {params.metadata.mAS} --Nt {params.metadata.Nt} --Ns {params.metadata.Ns}"
        " --min_trajectory {params.metadata.init_conf} --max_trajectory {params.metadata.final_conf} --trajectory_step {params.metadata.delta_conf_spectrum}"
        " --channel {wildcards.channel} --E0_plateau_start {params.plateau_start} --E0_plateau_end {params.plateau_end}"
 

def extraction_samples(wildcards):
    return [
        f"../JSONs/{dir_template}/meson_extraction_{rep}_{channel}_samples.json".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v", "t", "av", "at", "s"]
        for rep in ["f", "as"]
    ] + [
        f"../JSONs/{dir_template}/chimera_extraction_{channel}_{parity}_samples.json".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["lambda", "sigma", "sigmastar"]
        for parity in ["even", "odd"]
    ]

def mass_gevp_samples(wildcards):
    return [
        f"../JSONs/{dir_template}/meson_gevp_{rep}_{channel}_samples.json".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["ps", "v", "t", "av", "at", "s"]
        for rep in ["f", "as"]
    ] + [
        f"../JSONs/{dir_template}/chimera_gevp_{channel}_{parity}_samples.json".format(**row)
        for row in metadata.to_dict(orient="records")
        for channel in ["lambda", "sigma", "sigmastar"]
        for parity in ["even", "odd"]
    ]

rule mass_plot:
    params:
        module=lambda wildcards, input: input.script.replace("/", ".")[:-3],
    input:
        gevp_data=mass_gevp_samples,
        extraction_data=extraction_samples, #passing more data than the code needs to generate all the results. To be changed...
        script="src/plots/gevp_meson.py",
    output:
        plot="assets/plots/test.pdf",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python -m {params.module} {input.gevp_data} {input.extraction_data} --plot_file {output.plot}"
