using Pkg; Pkg.activate("."); Pkg.instantiate()
using DelimitedFiles
using MadrasSokal
println("Write table 1...")

files  = readdir("../../../../CSVs/output_topology/flow_analysis",join=true) 

ispath("../../../../CSVs/output_topology/") || mkpath("../../../../CSVs/output_topology/")
io1 = open("../../../../CSVs/output_topology/gradient_flow_observables.csv","w")
io3 = open("../../../../tables/table1.tex","w")

write(io1,"beta,T,L,mf,mas,w0,dw0,p,dp,Q,dQ,first,last,skip,Nconf\n")
#write(io2,"beta,mas,mf,Nt,Nl,first,skip,Nconf,p,w0,τ(Q),Q\n")

for file in files
    T, L, beta, mf,  mas = parse_filename(file)
    header  = readline(file)
    w0, dw0 = parse_w0(header)
    # get data for plaquette and topological charge
    data = readdlm(file,',',skipstart=1)
    traj = Int.(data[:,1])
    topo = data[:,2]
    plaq = data[:,3]
    # first and last configuration
    Nfirst, Nlast = extrema(traj)
    Nskip = traj[2] - traj[1]
    Nconf = length(traj)
    # get averages
    p, dp = stdmean(plaq;bin=2)
    Q, dQ = stdmean(topo;bin=2)
    # calculate autocorrelation times
    τP, dτP = madras_sokal_time(plaq)
    τQ, dτQ = madras_sokal_time(topo)
    # write to csv file
    #write(io1,"ensemble,beta,T,L,mf,mas,w0,dw0,p,dp,Q,dQ,first,last,skip,Nconf\n")
    #write(io2,"ensemble,beta,mas,mf,Nt,Nl,first,skip,Nconf,p,w0,tau(Q),Q\n")

    function ensemble_name(T,mf)
        T == 48 && return "M1"
        T == 64 && isapprox(mf,-0.71) && return "M2"
        T == 96 && return "M3"
        T == 64 && isapprox(mf,-0.70) && return "M4"
        T == 64 && isapprox(mf,-0.72) && return "M5"
    end
    

    # get ensemble name 
    MX = ensemble_name(T,mf)
    # write to csv file
    write(io1,"$MX,$beta,$T,$L,$mf,$mas,$w0,$dw0,$p,$dp,$Q,$dQ,$Nfirst,$Nlast,$Nskip,$Nconf\n")
    #write(io2,"$MX,$beta,$mas,$mf,$T,$L,$Nfirst,$Nskip,$Nconf,$(errorstring(p,Δp)),$(errorstring(ω0,Δω0)),$(errorstring(τQ,ΔτQ)),$(errorstring(Q,ΔQ))\n")
    write(io3,"$MX & $beta & $mas & $(rpad(mf,5,'0')) & $T & $L & $Nfirst & $Nskip & $Nconf & $(errorstring(p,dp)) & $(errorstring(w0,dw0)) & $(errorstring(τQ,dτQ)) & $(errorstring(Q,dQ))\\\\ \n")

    #write(io1,"$beta,$T,$L,$mf,$mas,$w0,$dw0,$p,$dp,$Q,$dQ,$Nfirst,$Nlast,$Nskip,$Nconf\n")
    #write(io3,"$beta & $mas & $mf & $T & $L & $Nfirst & $Nskip & $Nconf & $(errorstring(p,dp)) & $(errorstring(w0,dw0)) & $(errorstring(τQ,dτQ)) & $(errorstring(Q,dQ))\\\\ \n")
end

close(io1)
#close(io2)
close(io3)

