import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import os
from read_json import process_gevp_En_mass_samples
from argparse import ArgumentParser

def chimera_groundstate_plot(output_file):
    rep = 'as'
    ensembles = ['M1', 'M2', 'M3', 'M4', 'M5']
    json_dirs = [
        'JSONs/Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.71mAS-1.01T96L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.7mAS-1.01T64L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.72mAS-1.01T64L32',
    ]

    channels = ['Chimera_OC_even', 'Chimera_OV12_even']
    json_channels = ['lambda_even', 'sigma_even']
    tex_channel = ['$\\Lambda^{+}_{\\rm CB}$', '$\\Sigma^{+}_{\\rm CB}$']
    kernels = ['GAUSS','CAUCHY']
    ensembles_sep = 0.7
    kernel_spacing = 0.2
    peak_spacing = 0.1 

    stylefile = 'styles/paperdraft.mplstyle'
    plt.style.use(stylefile)

    colors = [cm.tab10(i) for i in np.linspace(0, 0.5, 6)]
    fig, ax = plt.subplots()

    for i in range(1,6):
        for o in range(2):
            for kernel in kernels:

                channel = channels[o]
                ensembles_offset = (i-1)*ensembles_sep

                CSV = f"CSVs/M{i}_chimerabaryons_spectral_density_spectrum.csv"
                df = pd.read_csv(CSV)
                df = df.loc[df['channel'] == channel]
                df = df.loc[df['rep'] == rep]
                data = df.loc[df['kernel'] == kernel]

                kernel_offset = 0 if kernel == 'GAUSS' else kernel_spacing
                peak_offset = peak_spacing * np.array([0, 1])
                shape = "o" if kernel == 'GAUSS' else "v"
                fs = 'full' if o == 1 else 'none'

                # get values to be plottted
                x = ensembles_offset + kernel_offset + peak_offset
                y = data['aE_0']
                Delta_y = data['errorE0']
            
                plt.errorbar(x, y, yerr=Delta_y, fmt=shape, color=colors[i%6], alpha=0.7, fillstyle=fs)
                plt.xlabel("ensemble")
                plt.ylabel("$a m_{CB}$")

                # add legend for marker shapes
                if i == 1:
                    plt.errorbar([], [], fmt=shape, color="k", alpha=0.7, label = f"{tex_channel[o]}: {kernel} kernel", fillstyle=fs)

            # get gevp data
            json_file = os.path.join(json_dirs[i-1],f"chimera_gevp_{json_channels[o]}_samples.json")
            E, DeltaE = process_gevp_En_mass_samples(json_file, json_channels[o], rep, 0)
            x_gevp = x[-1] + peak_spacing

            plt.errorbar([x_gevp], [E], yerr=DeltaE, fmt="D", color=colors[i%6], alpha=0.7, fillstyle=fs)
            if i == 1:
                plt.errorbar([], [], fmt="D" , color="k", alpha=0.7, label = f"{tex_channel[o]}: GEVP", fillstyle=fs)
            ax.legend(ncol=2)
    plt.xticks(ensembles_sep * np.array([0, 1, 2, 3, 4]) + kernel_spacing - peak_spacing/2, ['M1', 'M2', 'M3', 'M4', 'M5'])
    plt.savefig(output_file)


parser = ArgumentParser()
parser.add_argument("--outpdf",default=None, help="PDF file to be created")
args = parser.parse_args()
chimera_groundstate_plot(args.outpdf)
