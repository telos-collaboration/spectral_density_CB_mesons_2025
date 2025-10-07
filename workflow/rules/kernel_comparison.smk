import os


def input_files_new_meson_plot(rep,channel_label):
    d = [
        'JSONs/Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.71mAS-1.01T96L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.7mAS-1.01T64L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.72mAS-1.01T64L32',
    ]
    channel = f"f_{channel_label}" if rep == 'fund' else f"as_{channel_label}"
    json_files = [ os.path.join(d[i-1], f"meson_gevp_{channel}_samples.json") for i in range(1,6)]
    csv_files = [ f"CSVs/M{i}_spectral_density_spectrum.csv" for i in range(1,6) ]
    return json_files + csv_files


def input_files_new_chimera_plot():
    dirs = [
        'JSONs/Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.71mAS-1.01T96L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.7mAS-1.01T64L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.72mAS-1.01T64L32',
    ]
    states = ['lambda_even', 'sigma_even']

    return [ 
        os.path.join(dirs[i-1], f"chimera_gevp_{states[o]}_samples.json")
        for i in range(1,6) for o in range(2)        
    ]


rule chimera_kernel_plots:
    input:
        jsons=input_files_new_chimera_plot(),
        csvs=expand("CSVs/M{ensemble_idx}_chimerabaryons_spectral_density_spectrum.csv", ensemble_idx=[1, 2, 3, 4, 5]),
        script="just_plotting/code/kernel_comparison/plot_chimera.py",
    output:
        plot="assets/plots/kernel_comparison_chimera.pdf",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python {input.script} --outpdf {output.plot}" 


rule meson_kernel_plots:
    input:
        jsons=input_files_new_chimera_plot(),
        csvs=expand("CSVs/M{ensemble_idx}_spectral_density_spectrum.csv", ensemble_idx=[1, 2, 3, 4, 5]),
        script="just_plotting/code/kernel_comparison/plot_mesons.py",
    output:
        plot="assets/plots/kernel_comparison_mesons_{rep}_{channel}_{channel_label}.pdf",
    conda:
        "../envs/flow_analysis.yml"
    shell:
        "python {input.script} --outpdf {output.plot} --rep {wildcards.rep} --channel {wildcards.channel} --channel_label {wildcards.channel_label}"
