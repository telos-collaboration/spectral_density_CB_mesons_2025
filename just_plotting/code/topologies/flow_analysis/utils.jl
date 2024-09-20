using DelimitedFiles
function run_flow_analysis(logfile,outfile)
    try
        run(`python flow_analysis_cli.py $logfile $outfile`)
    catch
        run(`python3 flow_analysis_cli.py $logfile $outfile`)
    end
end
function plaquettes_log(file)
    plaquettes = Float64[]
    for line in eachline(file)
        if occursin("Plaquette",line)
            line = replace(line,"="=>" ")
            line = replace(line,":"=>" ")
            p = parse(Float64,split(line)[end])
            append!(plaquettes,p)
        end
    end
    return plaquettes
end
