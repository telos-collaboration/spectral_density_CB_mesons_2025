import lsdensities.utils.rhoUtils as u
from lsdensities.utils.rhoUtils import (
    init_precision,
    LogMessage,
    end,
    Inputs,
    create_out_paths,
    plot_markers,
    CB_colors,
    timesfont,
)
from lsdensities.utils.rhoParser import parseArgumentPrintSamples
from lsdensities.correlator.correlatorUtils import symmetrisePeriodicCorrelator
from lsdensities.utils.rhoParallelUtils import ParallelBootstrapLoop
from mpmath import mp, mpf
from lsdensities.core import A0E_mp, Smatrix_mp
from lsdensities.abw import gAg, gBg
from lsdensities.transform import h_Et_mp_Eslice, y_combine_sample_Eslice_mp_ToFile
import os
import time
from lsdensities.utils.rhoMath import invert_matrix_ge
import matplotlib.pyplot as plt
import csv
import read_hdf
import translate
#   Take input for Rho
import numpy as np

import os
import shutil

def copy_file(source_path, destination_path):
    shutil.copy(source_path, destination_path)

def create_out_paths2(path):
    if not os.path.exists(path):
        os.makedirs(path)

    plotpath = os.path.join(path, "Plots")
    logpath = os.path.join(path, "Logs")

    if not os.path.exists(plotpath):
        os.makedirs(plotpath)
    if not os.path.exists(logpath):
        os.makedirs(logpath)

    return plotpath, logpath

def init_variables_samples(datapath, outdir, ne, emin, emax, periodicity, kernel, sigma, prec, nboot, e0, Na, A0cut, mpi, rhopath):
    in_ = Inputs()
    in_.tmax = 0
    in_.periodicity = periodicity
    in_.kerneltype = kernel
    in_.prec = prec
    in_.datapath = datapath
    in_.outdir = outdir
    in_.massNorm = mpi
    in_.num_boot = nboot
    in_.sigma = sigma
    in_.emax = (
            emax * mpi
    )  # we pass it in unit of Mpi, here to turn it into lattice (working) units
    if emin == 0:
        in_.emin = (mpi / 20) * mpi
    else:
        in_.emin = emin * mpi
    in_.e0 = e0
    in_.Ne = ne
    in_.Na = Na
    in_.A0cut = A0cut
    in_.rhopath = rhopath
    return in_

def printSamples(datapath, outdir, ne, emin, emax, periodicity, kernel, sigma, prec, nboot, e0, Na, A0cut, mpi, rhopath, part_outdir):
    print(LogMessage(), "Initialising")
    #args = parseArgumentPrintSamples()
    init_precision(prec)
    par = init_variables_samples(datapath, outdir, ne, emin, emax, periodicity, kernel, sigma, prec, nboot, e0, Na, A0cut, mpi, rhopath)

    #   Reading datafile, storing correlator
    rawcorr, par.time_extent, par.num_samples = u.read_datafile(par.datapath)


    rho_file = rhopath
    inputrhofile = np.genfromtxt(rho_file, comments="#")
    energy = inputrhofile[:, 0]
    lambda_e = inputrhofile[:, 1]
    in_rho = inputrhofile[:, 2]
    in_stat = inputrhofile[:, 3]
    inputrhofile[:, 4]
    inputrhofile[:, 5]
    par.Ne = len(energy)
    espace = np.zeros(par.Ne)
    rho = np.zeros(par.Ne)
    drho = np.zeros(par.Ne)
    np.zeros(par.Ne)
    for _e in range(par.Ne):
        espace[_e] = energy[_e]

    #
    par.assign_values()
    par.report()
    par.plotpath, par.logpath = create_out_paths2(outdir)

    #   Here is the correlator
    rawcorr.evaluate()
    rawcorr.tmax = par.tmax

    #   Symmetrise
    if par.periodicity == "COSH":
        print(LogMessage(), "Folding correlator")
        symCorr = symmetrisePeriodicCorrelator(corr=rawcorr, par=par)
        symCorr.evaluate()

    #   Here is the resampling
    if par.periodicity == "EXP":
        corr = u.Obs(
            T=par.time_extent, tmax=par.tmax, nms=par.num_boot, is_resampled=True
        )
    if par.periodicity == "COSH":
        corr = u.Obs(
            T=symCorr.T,
            tmax=symCorr.tmax,
            nms=par.num_boot,
            is_resampled=True,
        )

    if par.periodicity == "COSH":
        # resample = ParallelBootstrapLoop(par, rawcorr.sample, is_folded=False)
        resample = ParallelBootstrapLoop(par, symCorr.sample, is_folded=False)
    if par.periodicity == "EXP":
        resample = ParallelBootstrapLoop(par, rawcorr.sample, is_folded=False)

    corr.sample = resample.run()
    corr.evaluate()

    #   Covariance
    print(LogMessage(), "Evaluate covariance")
    corr.evaluate_covmatrix(plot=False)
    corr.corrmat_from_covmat(plot=False)

    #   Make it into a mp sample
    print(LogMessage(), "Converting correlator into mpmath type")
    corr.fill_mp_sample()
    print(LogMessage(), "Cond[Cov C] = {:3.3e}".format(float(mp.cond(corr.mpcov))))
    cNorm = mpf(str(corr.central[1] ** 2))

    # from HLT class
    A0set = A0E_mp(espace, par, alpha_=0, e0_=par.mpe0)

    for _e in range(par.Ne):
        estar_ = espace[_e]
        fname = "lsdensitiesamplesE" + str(estar_) + "sig" + str(par.sigma)
        fpath = os.path.join(par.logpath, fname)
        _Bnorm = cNorm / (estar_ * estar_)
        _factor = (lambda_e[_e] * A0set[_e]) / _Bnorm
        S_ = Smatrix_mp(
            tmax_=par.tmax,
            alpha_=0,
            e0_=par.mpe0,
            type=par.periodicity,
            T=par.time_extent,
        )
        _M = S_ + (_factor * corr.mpcov)
        start_time = time.time()
        _Minv = invert_matrix_ge(_M)
        end_time = time.time()
        print(
            LogMessage(),
            "\t \t lambdaToRho ::: Matrix inverted in {:4.4f}".format(
                end_time - start_time
            ),
            "s",
        )
        _g_t_estar = h_Et_mp_Eslice(_Minv, par, estar_, alpha_=0)
        rho[_e], drho[_e] = y_combine_sample_Eslice_mp_ToFile(
            fpath, _g_t_estar, corr.mpsample, par
        )

        gag_estar = gAg(S_, _g_t_estar, estar_, 0, par)

        gBg_estar = gBg(_g_t_estar, corr.mpcov, _Bnorm)

        print(LogMessage(), "\t \t  B / Bnorm = ", float(gBg_estar))
        print(LogMessage(), "\t \t  A / A0 = ", float(gag_estar / A0set[_e]))

    plt.errorbar(
        x=espace,
        y=rho,
        yerr=drho,
        marker=plot_markers[0],
        markersize=2,
        elinewidth=1.3,
        capsize=2,
        ls="",
        label=r"Recomputed",
        color=CB_colors[0],
    )
    plt.errorbar(
        x=espace,
        y=in_rho,
        yerr=in_stat,
        marker=plot_markers[1],
        markersize=2,
        elinewidth=1.3,
        capsize=2,
        ls="",
        label=r"Input",
        color=CB_colors[1],
    )
    plt.title(r" $\sigma$" + " = {:2.2f} Mpi".format(par.sigma / par.massNorm))
    plt.xlabel(r"$E / M_{\pi}$", fontdict=timesfont)
    # plt.ylabel("Spectral density", fontdict=u.timesfont)
    plt.legend(prop={"size": 12, "family": "Helvetica"})
    plt.grid()
    plt.tight_layout()
    source_path = rhopath
    destination_path = part_outdir
    copy_file(source_path, destination_path)
    #plt.show()
    #end()

def main():
    file_path = '../input_correlators/chimera_data_full.hdf5'
    ####################### External data for make rho finding easier #######################
    categories = ['PS', 'V', 'T', 'AV', 'AT', 'S', 'ps', 'v', 't', 'av', 'at', 's']
    # Mesonic channels
    mesonic_channels = ['Chimera_OC_even', 'Chimera_OC_odd', 'Chimera_OV12_even', 'Chimera_OV12_odd', 'Chimera_OV32_even', 'Chimera_OV32_odd']
    # Ensembles: M1, M2, M3, M4, M5
    ensembles = ['M1', 'M2', 'M3', 'M4', 'M5']
    #ensembles = ['M1']
    # Roots in HDF5 for each ensemble
    roots = ['chimera_out_48x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1',
             'chimera_out_64x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1',
             'chimera_out_96x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1',
             'chimera_out_64x20x20x20nc4nf2nas3b6.5mf0.70mas1.01_APE0.4N50_smf0.2as0.12_s1',
             'chimera_out_64x32x32x32nc4nf2nas3b6.5mf0.72mas1.01_APE0.4N50_smf0.24as0.12_s1']
    # Representations considered
    #rep = 'fund'
    reps = ['fund', 'anti']
    # Kernel in HLT
    kerneltype = ['HALFNORMGAUSS', 'CAUCHY']
    #kerneltype = ['HALFNORMGAUSS']
    # Initialize dictionaries to store the data
    Nsource_C_values_MN = {}
    Nsink_C_values_MN = {}
    am_C_values_MN = {}
    sigma1_over_mC_values_MN = {}
    sigma2_over_mC_values_MN = {}
    # Read data from CSV
    with open('metadata/metadata_spectralDensity_chimerabaryons.csv', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ensemble = row['Ensemble']
            # Initialize lists for each ensemble if not already present
            if ensemble not in Nsource_C_values_MN:
                Nsource_C_values_MN[ensemble] = []
                Nsink_C_values_MN[ensemble] = []
                am_C_values_MN[ensemble] = []
                sigma1_over_mC_values_MN[ensemble] = []
                sigma2_over_mC_values_MN[ensemble] = []
            # Append data for each category to the respective lists
            for category in categories:
                Nsource_C_values_MN[ensemble].append(int(row[f'{category}_Nsource']))
                Nsink_C_values_MN[ensemble].append(int(row[f'{category}_Nsink']))
                am_C_values_MN[ensemble].append(float(row[f'{category}_am']))
                sigma1_over_mC_values_MN[ensemble].append(float(row[f'{category}_sigma1_over_m']))
                sigma2_over_mC_values_MN[ensemble].append(float(row[f'{category}_sigma2_over_m']))
    # Create a 3D matrix with ensemble index
    matrix_4D = [
        [
            ensemble,
            am_C_values_MN[ensemble],
            sigma1_over_mC_values_MN[ensemble],
            sigma2_over_mC_values_MN[ensemble],
            Nsource_C_values_MN[ensemble],
            Nsink_C_values_MN[ensemble]
        ]
        for ensemble in ensembles
    ]

    for rep in reps:
        for kernel in kerneltype:
            for index, ensemble in enumerate(ensembles):
                for k, channel in enumerate(mesonic_channels):
                    if rep == 'fund':
                        Nsource = matrix_4D[index][4][k]
                        Nsink = matrix_4D[index][5][k]
                    else:
                        Nsource = matrix_4D[index][4][k + 6]
                        Nsink = matrix_4D[index][5][k + 6]
                    group_prefixes = {
                        'gi': ['g1', 'g2', 'g3'],
                        'g0gi': ['g0g1', 'g0g2', 'g0g3'],
                        'g5gi': ['g5g1', 'g5g2', 'g5g3'],
                        'g0g5gi': ['g0g5g1', 'g0g5g2', 'g0g5g3']
                    }
                    prefix = group_prefixes.get(channel, [channel])
                    datasets = []
                    for g in prefix:
                        dataset_path = f"{roots[index]}/source_N{Nsource}_sink_N{Nsink}/{channel}_re"
                        group1 = f"source_N{Nsource}_sink_N{Nsink}"
                        group2 = f"{channel}_re"
                        datasets.append(read_hdf.extract_dataset(file_path, group2, roots[index], group1))
                        with open('paths.log', 'a') as file:
                            print(dataset_path, file=file)
                    dataset = sum(datasets) / len(datasets) if datasets else None
                    if dataset is not None:
                        translate.save_matrix_to_file(dataset, f'corr_to_analyse_{channel}_{rep}_{ensemble}.txt')
                    mpi = matrix_4D[index][1][k]
                    if kernel == 'HALFNORMGAUSS':
                        tmp = mpi * matrix_4D[index][2][k]
                    elif kernel == 'CAUCHY':
                        tmp = mpi * matrix_4D[index][3][k]
                    mpi = mpi
                    sigma = tmp
                    decimal_part = tmp / matrix_4D[index][1][k] % 1
                    decimal_as_int = int(decimal_part * 100)
                    datapath = f'./corr_to_analyse_{channel}_{rep}_{ensemble}.txt'
                    if kernel == 'HALFNORMGAUSS':
                        kernel2 = 'GAUSS'
                    elif kernel == 'CAUCHY':
                        kernel2 = kernel
                    spdens_outdir = f'./{ensemble}_{channel}_s0p{decimal_as_int}_{kernel}_Nsource{Nsource}_Nsink{Nsink}'
                    all_items = os.listdir(spdens_outdir)
                    subdirectories = [item for item in all_items if os.path.isdir(os.path.join(spdens_outdir, item))]
                    if len(subdirectories) == 1:
                        subdirectory_name = subdirectories[0]
                        spdens_outdir = os.path.join(spdens_outdir, subdirectory_name)
                        print("Path to the only subdirectory:", spdens_outdir)
                    else:
                        print("There is either no subdirectory or more than one subdirectory in the directory.")
                    spdens_outdir = spdens_outdir + '/Logs/ResultHLT.txt'

                    
                    outdir = f'../input_fit/{ensemble}/{channel}_Nsource{Nsource}_Nsink{Nsink}/{kernel2}/{channel}_Nsource{Nsource}_Nsink{Nsink}'
                    part_outdir = f'../input_fit/{ensemble}/{channel}_Nsource{Nsource}_Nsink{Nsink}/{kernel2}/fit_results.txt'

                    ne = 15
                    emin = 0.3
                    emax = 2.8
                    periodicity = 'COSH'
                    prec = 105
                    nboot = 300
                    e0 = 0.0
                    Na = 1
                    A0cut = 0.1


                    printSamples(datapath, outdir, ne, emin, emax, periodicity, kernel, sigma, prec, nboot, e0, Na, A0cut, mpi, spdens_outdir, part_outdir)


if __name__ == "__main__":
    main()
