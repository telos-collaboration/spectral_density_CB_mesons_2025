from flow_analysis.readers.read_hirep import read_flows_hirep
from flow_analysis.measurements.scales import measure_w0
import sys

def flow_hirep(hirep_file, W0, filename):
    # obtain topological charge Q, the number of the configurations as specified in the configuration filenames and the
    # gradient flow scale. For the ensembles in this dataset the reference scale has been fixed to W0=0.28125 in the
    # measurements.
    flows = read_flows_hirep(hirep_file)
    Qs = flows.Q_history()
    trajectories = flows.trajectories
    scale = measure_w0(flows, W0)

    # Write everything into a csv with that is compatible with the current julia scripts used in the analysis
    f = open(filename, "w")
    f.write("trajectory,Q (w0 = %s)\n" % scale)
    for i in zip(trajectories, Qs):
        f.write("%s,%s\n" % (i[0], i[1]))
    f.close()

args = sys.argv
if len(args) < 3:
    print("Missing input and/or output file")
elif len(args)==3:
    hirep_file  = args[1]
    output_file = args[2] 
    W0 = 0.28125
    flow_hirep(hirep_file, W0, output_file)
elif len(args)>3:
    hirep_file  = args[1]
    output_file = args[2] 
    W0 = args[3]
    flow_hirep(hirep_file, W0, output_file)
