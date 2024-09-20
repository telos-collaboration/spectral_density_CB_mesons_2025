function plaquettes_grid(file)
    plaquettes = Float64[]
    configurations = Int64[]
    for line in eachline(file)
        p = findfirst("Plaquette",line)
        line = line[p[end]+2:end]
        line = replace(line,"[" =>" ")
        line = replace(line,"]" =>" ")
        data =  split(line)
        append!(configurations,parse(Int64,data[1]))
        append!(plaquettes,parse(Float64,data[2]))
    end
    perm = sortperm(configurations)
    permute!(plaquettes,perm)
    permute!(configurations,perm)
    # only keep one value for every configurations
    # unique indices
    _unique_indices(x) = unique(i -> x[i],1:length(x))
    inds = _unique_indices(configurations)
    configurations = getindex(configurations,inds)
    plaquettes = getindex(plaquettes,inds)
    return configurations, plaquettes
end
function _parameters_from_filename(file)
    info = split(replace(file,r"[A-Za-z,_,//]+" => " "))[2:end]
    T, L = parse.(Int,info[1:2])
    β, amf0, amas0 = parse.(Float64,info[3:end]) 
    amf0  *= -1
    amas0 *= -1
    return T, L, β, amf0, amas0
end
function _wilson_flow_scale_from_file(file)
    header    = readline(file)
    w0_string = replace(header,"=" => " ", "+/-" => " ", r"[(,)]"=>" ")
    w0, dw0   = parse.(Float64,split(w0_string)[end-1:end])
    return w0, dw0
end
function _binned_mean_std(O;bin=1)
    mO = mean(O)
    sO = std(O)/sqrt( length(O) / bin )
    return mO, sO
end 
