using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("utils.jl")

inputDIR  = abspath("../../../../input_topology/gradient_flow/")
outputDIR = abspath("../../../../CSVs/output_topology/flow_analysis/")
ispath(outputDIR) || mkpath(outputDIR)
ensembles = readdir(inputDIR)

println("Analyse gradient flow output files")
for ensemble in ensembles
    @show ensemble
 
    hirep_file = joinpath(inputDIR,"$ensemble/out_flow")
    output_file = joinpath(outputDIR,ensemble*"_flow")
    run_flow_analysis(hirep_file,output_file)

    # now get the plaquette
    plaq = plaquettes_log(hirep_file)
    data = readdlm(output_file,',',skipstart=1)
    head = readline(output_file)

    io = open(output_file,"w")
    write(io,head*",plaquette\n")
    writedlm(io,hcat(data,plaq),',')
    close(io)
end
