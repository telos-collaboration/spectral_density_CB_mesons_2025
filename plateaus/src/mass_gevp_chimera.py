#!/usr/bin/env python3


import h5py
import numpy as np
import logging

from .bootstrap import BootstrapSampleSet, bootstrap_finalize, BOOTSTRAP_SAMPLE_COUNT
from .dump import dump_dict, dump_samples
from . import extract
from .mass import (
    get_baryon_corr,
    get_args,
)
from .read_hdf5 import get_ensemble
from .plots_common import plot_baryon_gevp_energy_states


def get_baryon_Cmat_single(ensemble, args, Nmin, Nmax, Nd, channel):
    source_Ns = [str(i) for i in np.arange(Nmin, Nmax + 1, Nd)]
    sink_Ns = [str(i) for i in np.arange(Nmin, Nmax + 1, Nd)]

    size = len(source_Ns)
    

    mat_e_sample = np.zeros(shape=(BOOTSTRAP_SAMPLE_COUNT, args.Nt, size, size))
    mat_o_sample = np.zeros(shape=(BOOTSTRAP_SAMPLE_COUNT, args.Nt, size, size))
    mat_e_mean = np.zeros(shape=(1, args.Nt, size, size))
    mat_o_mean = np.zeros(shape=(1, args.Nt, size, size))

    # print(mat.shape)

    for i in range(size):
        for j in range(size):

            corr_e, corr_o = get_baryon_corr(
                ensemble, args, source_Ns[i], sink_Ns[j], channel
                )
            
            mat_e_sample[:, :, i, j] = corr_e.samples
            mat_o_sample[:, :, i, j] = corr_o.samples
            mat_e_mean[:, :, i, j] = corr_e.mean
            mat_o_mean[:, :, i, j] = corr_o.mean

    return BootstrapSampleSet(mat_e_mean, mat_e_sample), BootstrapSampleSet(mat_o_mean, mat_o_sample)


def gevp_extraction(ensemble, args):

    corr_mat_even, corr_mat_odd = get_baryon_Cmat_single(ensemble, args, 0, 80, 40, args.channel)
    if args.channel.split("_")[1] == "even":
        corr_mat = corr_mat_even
    elif args.channel.split("_")[1] == "odd":
        corr_mat = corr_mat_odd
    else:
        logging.waring("NO such parity: ")

    eigenvalues = extract.gevp_fixT(corr_mat.mean, corr_mat.samples, args.gevp_t0, args.gevp_t0+1, args.Nt)

    return eigenvalues


def main():
    args = get_args()

    data = h5py.File(args.h5file, "r")
    ensemble, = get_ensemble(
        data,
        beta=args.beta,
        mF=args.mF,
        mAS=args.mAS,
        Nt=args.Nt,
        Ns=args.Ns,
        epsilon=args.epsilon,
    )

    eigenvalues = gevp_extraction(ensemble, args)
    masses, chiquares = extract.extract_energy_states(eigenvalues, args)

    if args.effmass_plot_file:
        plot_baryon_gevp_energy_states(args, eigenvalues, masses)

    metadata = {
        "ensemble_name": args.ensemble_name,
        "beta": args.beta,
        "mAS": args.mAS,
        "Nt": args.Nt,
        "Ns": args.Ns,
    }

    if args.output_file_samples:
        data_to_save = {**metadata}

        for n, mass in enumerate(masses):
            data_to_save[f"gevp_{args.channel}_E{n}_chisquare"] = chiquares[n]
            data_to_save[f"gevp_{args.channel}_E{n}_mass"] = mass

        dump_samples(data_to_save,args.output_file_samples)


if __name__ == "__main__":
    main()
