import numpy as np
import h5py
import csv
import sys
import matplotlib.pyplot as plt
import os.path
import pandas as pd

sys.path.insert(1, "./Lib")

import extract
import bootstrap
import read_hdf
import fitting


ens_lb = {
    "M1": "chimera_out_48x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "M2": "chimera_out_64x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "M3": "chimera_out_96x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "M4": "chimera_out_64x20x20x20nc4nf2nas3b6.5mf0.70mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "M5": "chimera_out_64x32x32x32nc4nf2nas3b6.5mf0.72mas1.01_APE0.4N50_smf0.24as0.12_s1",
}

DATA = h5py.File("../input_correlators/chimera_data_full.hdf5")

CSV_data = []
CHs = [
    ["g5"],
    ["g1", "g2", "g3"],
    ["g0g1", "g0g2", "g0g3"],
    ["g5g1", "g5g2", "g5g3"],
    ["g0g5g1", "g0g5g2", "g0g5g3"],
    ["id"],
]

CHs_tag = ["g5", "gi", "g0gi", "g5gi", "g0g5gi", "id"]


REP = "anti"

ENS = [ens_lb["M3"]]

CHs = [
    ["g1", "g0g1"],
    ["g2", "g0g2"],
    ["g3", "g0g3"],
]

En = 0

for ens_tag in ENS:
    tmp = ens_tag.split("_")[2]
    ens = tmp.split("nc")[0] + "b" + tmp.split("b")[1]

    print(ens)

    ens_group = DATA[ens_tag]

    Nt = int(ens.split("x")[0])
    Ns = int(ens.split("x")[1])
    beta = float(ens.split("b")[1].split("m")[0])
    as_bare_mass = ens_group["quarkmasses_antisymmetric"][0]
    epsilon_as = ens_group["Wuppertal_eps_anti"][0]

    t0_GEVP = 1  # t0 for GEVP
    print("\n t0 for GEVP =", t0_GEVP)

    tmp_bin = []

    for i in range(len(CHs)):
        ch = CHs[i]
        tmp_bin.append(
            read_hdf.get_meson_Cmat_mix(DATA, ens_tag, REP, 0, 80, 40, ch[0], ch[1])
        )

    Cmat = np.array(tmp_bin).mean(axis=0)

    LAM, VEC = extract.GEVP_fixT(
        Cmat, t0_GEVP, t0_GEVP + 1, Nt
    )  # t0, t1, t2 --> t0 for the T_fix in GEVP. t1 and t2 are the range of sloving GEVP.

    M_tmp = extract.Analysis_Mass_eff_simple(LAM[:, :, En], 1, Nt, 1, f"E{En}")

    plt.ylim(0.3, 2)
    plt.xlim(0, Nt / 2)

    plt.legend()
    plt.show()

    ti = int(input(f"E{En}: ti = "))
    tf = int(input(f"E{En}: tf = "))

    for tt in np.arange(ti + 3, tf + 1):
        m_tmp, m_tmp_err, chi2 = fitting.fit_exp_std(
            LAM[:, :, En],
            ti,
            tt,
            Nt,
        )
