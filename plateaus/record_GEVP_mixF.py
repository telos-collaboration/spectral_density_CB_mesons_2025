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


ens_lb = {
    "chimera_out_48x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1": "M1",
    "chimera_out_64x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1": "M2",
    "chimera_out_96x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1": "M3",
    "chimera_out_64x20x20x20nc4nf2nas3b6.5mf0.70mas1.01_APE0.4N50_smf0.2as0.12_s1": "M4",
    "chimera_out_64x32x32x32nc4nf2nas3b6.5mf0.72mas1.01_APE0.4N50_smf0.24as0.12_s1": "M5",
}

DATA = h5py.File("../input_correlators/chimera_data_full.hdf5")

CSV_data = []
CHs = [
    ["g1", "g0g1"],
    ["g2", "g0g2"],
    ["g3", "g0g3"],
]

savename = "metadata/MIX_GEVP_meta_F.csv"

n = 0
if os.path.isfile(savename):
    meta = pd.read_csv(savename)
    existence = True

else:
    meta = []
    existence = False

log_name = DATA.keys()


N_patience = 1

for ens_tag in log_name:
    tmp = ens_tag.split("_")[2]
    ens = tmp.split("nc")[0] + "b" + tmp.split("b")[1]

    if existence and ens in meta.ENS.values:
        print(ens, "already done")
        continue

    your_patience = int(input("How many ensemble you want to check? "))

    print(ens)

    ens_group = DATA[ens_tag]

    Nt = int(ens.split("x")[0])
    Ns = int(ens.split("x")[1])
    beta = float(ens.split("b")[1].split("m")[0])
    f_bare_mass = ens_group["quarkmasses_fundamental"][0]
    epsilon_f = ens_group["Wuppertal_eps_fund"][0]

    t0_GEVP = 1  # t0 for GEVP
    print("\n t0 for GEVP =", t0_GEVP)

    data_tmp = []
    data_tmp.extend((ens, Nt, Ns, beta, f_bare_mass, epsilon_f, t0_GEVP))

    tmp_bin = []
    for i in range(len(CHs)):
        ch = CHs[i]

        tmp_bin.append(
            read_hdf.get_meson_Cmat_mix(DATA, ens_tag, "fund", 0, 80, 40, ch[0], ch[1])
        )

    Cmat = np.array(tmp_bin).mean(axis=0)

    LAM, VEC = extract.GEVP_fixT(
        Cmat, t0_GEVP, t0_GEVP + 1, Nt
    )  # t0, t1, t2 --> t0 for the T_fix in GEVP. t1 and t2 are the range of sloving GEVP.

    for n in range(Cmat.shape[-1]):
        M_tmp = extract.Analysis_Mass_eff_simple(LAM[:, :, n], 1, Nt - 1, 1, f"E{n}")

        plt.ylim(0.3, 2)
        plt.xlim(0, Nt / 2)
        extract.sperater.reset()

        plt.legend()

        plt.show()

        ti = int(input(f"E{n}: ti = "))
        tf = int(input(f"E{n}: tf = "))
        data_tmp.extend((ti, tf))

    CSV_data.append(data_tmp)

    if N_patience + 1 > your_patience:
        break
    N_patience += 1


if existence:
    with open(savename, "a", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(CSV_data)

else:
    with open(savename, "w", newline="") as csvfile:
        fieldnames = [
            "ENS",
            "Nt",
            "Ns",
            "beta",
            "f_bare_mass",
            "f_epsilon",
            "t0_GEVP",
            "E0_ti",
            "E0_tf",
            "E1_ti",
            "E1_tf",
            "E2_ti",
            "E2_tf",
            "E3_ti",
            "E3_tf",
            "E4_ti",
            "E4_tf",
            "E5_ti",
            "E5_tf",
        ]

        writer = csv.writer(csvfile)
        writer.writerow(fieldnames)
        writer.writerows(CSV_data)
