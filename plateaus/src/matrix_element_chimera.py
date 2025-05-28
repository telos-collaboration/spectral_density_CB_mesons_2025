#!/usr/bin/env python3


import h5py
import numpy as np
import logging

from .bootstrap import BootstrapSampleSet, bootstrap_finalize, BOOTSTRAP_SAMPLE_COUNT
from .dump import dump_dict, dump_samples
from . import extract
from .fitting import fit_exp_simultaneous
from .mass import (
    get_baryon_corr,
    get_args,
)
from .read_hdf5 import get_ensemble
from .plots_common import plot_baryon_gevp_energy_states




def baryon_extraction(ensemble, args):

    corr_e_ss, corr_o_ss = get_baryon_corr(ensemble, args, 80, 80, args.channel)
    corr_e_sp, corr_o_sp = get_baryon_corr(ensemble, args, 80, 0, args.channel)

    parity = args.channel.split("_")[1]
    if parity == "even":
        corr_ss, corr_sp = corr_e_ss, corr_e_sp
    elif parity == "odd":
        corr_ss, corr_sp = corr_o_ss, corr_o_sp
    else:
        logging.waring(f"NO such parity: {parity}")

    mass, matrix_element, chi2 = fit_exp_simultaneous(
        corr_ss,
        corr_sp,
        args.E0_plateau_start,
        args.E0_plateau_end,
    )
    result_samples= matrix_element.samples *np.sqrt(2*mass.samples)
    return mass, result_samples, chi2


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
    
    mass, matrix_element, chisquare = baryon_extraction(ensemble, args)
    
    metadata = {
        "ensemble_name": args.ensemble_name,
        "beta": args.beta,
        "mAS": args.mAS,
        "Nt": args.Nt,
        "Ns": args.Ns,
    }

    if args.output_file_samples:
        data_to_save = {**metadata}

        data_to_save[f"{args.channel}_chisquare"] = chisquare
        data_to_save[f"{args.channel}_mass"] = mass
        data_to_save[f"{args.channel}_matrix_element"] = matrix_element

        dump_samples(data_to_save,args.output_file_samples)


if __name__ == "__main__":
    main()
