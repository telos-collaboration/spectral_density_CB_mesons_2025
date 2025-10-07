rule parse_chimera_smeared:
    input:
        files=glob("raw_data/smear/*"),
        script="plateaus/write_chimera.jl",
        metadata="metadata/channels_chimera.txt"
    output:
        h5="input_correlators/chimera_data_reduced.h5",
    conda:
        "../envs/hirep_parsing.yml"
    shell:
        'julia --project="libs/HiRepParsing" {input.script} --h5file {output.h5} --channels {input.metadata} {input.files}'
