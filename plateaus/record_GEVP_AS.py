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
    ["g5"],
    ["g1", "g2", "g3"],
    ["g0g1", "g0g2", "g0g3"],
    ["g5g1", "g5g2", "g5g3"],
    ["g0g5g1", "g0g5g2", "g0g5g3"],
    ["id"],
]

CHs_tag = ["g5", "gi", "g0gi", "g5gi", "g0g5gi", "id"]


savename = "metadata/mesonAS_GEVP_meta.csv"

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

    your_patience = int(input("How many ensemble you want to check?"))

    print(ens)

    ens_group = DATA[ens_tag]

    Nt = int(ens.split("x")[0])
    Ns = int(ens.split("x")[1])
    beta = float(ens.split("b")[1].split("m")[0])
    as_bare_mass = ens_group["quarkmasses_antisymmetric"][0]
    epsilon_as = ens_group["Wuppertal_eps_anti"][0]

    t0_GEVP = 1  # t0 for GEVP
    print("\n t0 for GEVP =", t0_GEVP)

    data_tmp = []
    data_tmp.extend((ens, Nt, Ns, beta, as_bare_mass, epsilon_as, t0_GEVP))

    for i in range(len(CHs)):
        ch = CHs[i]

        tmp_bin = []

        print(CHs_tag[i])

        for j in range(len(ch)):
            tmp_bin.append(
                read_hdf.get_meson_Cmat_single(DATA, ens_tag, "anti", 0, 80, 40, ch[j])
            )

        Cmat = np.array(tmp_bin).mean(axis=0)

        LAM, VEC = extract.GEVP_fixT(
            Cmat, t0_GEVP, t0_GEVP + 1, Nt
        )  # t0, t1, t2 --> t0 for the T_fix in GEVP. t1 and t2 are the range of sloving GEVP.

        for n in range(Cmat.shape[-1]):
            M_tmp = extract.Analysis_Mass_eff_cosh(LAM[:, :, n], 1, Nt, f"E{n}")

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
            "as_bare_mass",
            "as_epsilon",
            "t0_GEVP",
            "g5_E0_ti",
            "g5_E0_tf",
            "g5_E1_ti",
            "g5_E1_tf",
            "g5_E2_ti",
            "g5_E2_tf",
            "gi_E0_ti",
            "gi_E0_tf",
            "gi_E1_ti",
            "gi_E1_tf",
            "gi_E2_ti",
            "gi_E2_tf",
            "g0gi_E0_ti",
            "g0gi_E0_tf",
            "g0gi_E1_ti",
            "g0gi_E1_tf",
            "g0gi_E2_ti",
            "g0gi_E2_tf",
            "g5gi_E0_ti",
            "g5gi_E0_tf",
            "g5gi_E1_ti",
            "g5gi_E1_tf",
            "g5gi_E2_ti",
            "g5gi_E2_tf",
            "g0g5gi_E0_ti",
            "g0g5gi_E0_tf",
            "g0g5gi_E1_ti",
            "g0g5gi_E1_tf",
            "g0g5gi_E2_ti",
            "g0g5gi_E2_tf",
            "id_E0_ti",
            "id_E0_tf",
            "id_E1_ti",
            "id_E1_tf",
            "id_E2_ti",
            "id_E2_tf",
        ]

        writer = csv.writer(csvfile)
        writer.writerow(fieldnames)
        writer.writerows(CSV_data)
