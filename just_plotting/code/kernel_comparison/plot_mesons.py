import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import os
from read_json import process_gevp_En_mass_samples
from argparse import ArgumentParser

def comparison_plot_kernels(output_file,rep,channel,channel_label):

    stylefile = 'styles/paperdraft.mplstyle'
    kernels = ['GAUSS','CAUCHY']
    json_channel = f"f_{channel_label}" if rep == 'fund' else f"as_{channel_label}"

    ensembles_sep = 0.7
    kernel_spacing = 0.2
    peak_spacing = 0.1 

    colors = [cm.tab10(i) for i in np.linspace(0, 0.5, 6)]
    fig, ax = plt.subplots()

    ensembles = ['M1', 'M2', 'M3', 'M4', 'M5']
    json_dirs = [
        'JSONs/Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.71mAS-1.01T96L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.7mAS-1.01T64L20',
        'JSONs/Sp4b6.5nF2nAS3mF-0.72mAS-1.01T64L32',
        ]

    for i in range(1,6):
        for kernel in kernels:

            ensembles_offset = (i-1)*ensembles_sep

            CSV = f"CSVs/M{i}_spectral_density_spectrum.csv"
            df = pd.read_csv(CSV)
            df = df.loc[df['channel'] == channel]
            df = df.loc[df['rep'] == rep]
            data = df.loc[df['kernel'] == kernel]

            kernel_offset = 0 if kernel == 'GAUSS' else kernel_spacing
            peak_offset = peak_spacing * np.array([0, 1])
            shape = "o" if kernel == 'GAUSS' else "v"

            # get values to be plottted
            x = ensembles_offset + kernel_offset + peak_offset
            y = np.array(data['aE_0'])
            Delta_y = np.array(data['errorE0'])

            # Sometimes there are no uncertainties given by lsd
            # Niccolò sets the error to be 1%. Let me add 1% of 
            # systematics to the uncertainties.
            for j,val in enumerate(Delta_y):
                if Delta_y[j] < 1E-10:
                   Delta_y[j] = y[j]/100 
        
            ylabel = '$a m_{\\rm PS}$' if rep == "fund" else '$a m_{\\rm ps}$'
            plt.errorbar(x, y, yerr=Delta_y, fmt=shape, color=colors[i%6], alpha=0.7)
            plt.ylabel(ylabel)
            plt.xlabel("ensemble")

            # add legend for marker shapes
            if i == 1:
                plt.errorbar([], [], fmt=shape, color="k", alpha=0.7, label = f"{kernel} kernel")

        # get gevp data
        json_file = os.path.join(json_dirs[i-1],f"meson_gevp_{json_channel}_samples.json")
        E, DeltaE = process_gevp_En_mass_samples(json_file, json_channel, rep, 0)
        x_gevp = x[-1] + peak_spacing

        plt.errorbar([x_gevp], [E], yerr=DeltaE, fmt="D", color=colors[i%6], alpha=0.7)
        if i == 1:
            plt.errorbar([], [], fmt="D" , color="k", alpha=0.7, label = "GEVP")
        ax.legend()

    # Plotting parameters
    plt.xticks(ensembles_sep * np.array([0, 1, 2, 3, 4]) + kernel_spacing - peak_spacing/2, ['M1', 'M2', 'M3', 'M4', 'M5'])
    plt.style.use(stylefile)
    plt.savefig(output_file)


parser = ArgumentParser()
parser.add_argument("--rep",default=None, help="Fermion representation to be analysed")
parser.add_argument("--channel",default=None, help="Meson channel to be analysed")
parser.add_argument("--channel_label",default=None, help="Label of meson channel according to definitions in the paper")
parser.add_argument("--outpdf",default=None, help="PDF file to be created")
args = parser.parse_args()

comparison_plot_kernels(args.outpdf,args.rep,args.channel,args.channel_label)