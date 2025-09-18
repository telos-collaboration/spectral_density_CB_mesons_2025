using Pkg
using HiRepParsing
using HDF5
using ArgParse

# This script parses the log files in the directory 'dir', and saves them as an hdf5-file 
# in the location provided by 'h5file'.

# It creates a single hdf5 file for all log files. Measurements performed on the same ensemble
# are written in distinct hdf5 groups labelled  by the variable `ensemble`

function main(files,h5file;filter_channels=false,channels=nothing,setup=true)

    # I have defined a reference variable for the last saved ensemble.
    # If the ensemble changes, we save also information on the lattice setup (coupling, size, bare masses) 
    # to the hdf5 file. (This is not very robust and depends o the naming scheme of the output files)
    ensemble0 = "" 
    for file in sort(files)

        # set up a regular expression, that matches the measurement type and the different smearing levels.
        # Note that the steps in the smearing levels are hard-coded
        regex = r"N(?<N1>[0-9]+)_N(?<N2>[0-9]+)"
        m = match(regex,basename(file))
        N1 = parse(Int,m[:N1])
        N2 = parse(Int,m[:N2])
        name  = replace(basename(file),r"(\.txt|\.zst|\.zstd)"=>"")
        
        # parse the ensemble name from the filename 
        # (again this depends strongly on the naming scheme)
        # Check if we need to write the lattice-setupparameters to hdf5 file
        ensemble = replace(name,"N$(N1)_N$(N2)"=>"")
        if !isempty(ensemble0)
            setup = ensemble != ensemble0
        end
        ensemble0 = ensemble
              
        regex = r"source_N[0-9]+_sink_N(0|40|80)"
        writehdf5_spectrum_with_regexp(file,h5file,regex;mixed_rep=true,h5group=ensemble,setup,filter_channels,channels,sort=true)
    end
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the parsed data"
        required = true
        "--channels"
        help = "File containing the channels to be parsed"
        required = true
        "files"
        nargs = '+'
        help = "Files to be parsed."
        required = true
    end
    return parse_args(s)
end

args = parse_commandline()
files = args["files"]
h5file = args["h5file"]
channels = readlines(args["channels"])
filter_channels=true

main(files,h5file;filter_channels,channels)