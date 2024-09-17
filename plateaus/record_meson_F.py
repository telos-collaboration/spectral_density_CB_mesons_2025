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


def LB_chimera(b):
    if b == "Lambda":
        return "Chimera_OC"
    elif b == "Sigma":
        return "Chimera_OV12"
    elif b == "SigmaS":
        return "Chimera_OV32"
    else:
        return "Chimera"


ens_lb = {
    "chimera_out_48x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1": "M1",
    "chimera_out_64x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1": "M2",
    "chimera_out_96x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1": "M3",
    "chimera_out_64x20x20x20nc4nf2nas3b6.5mf0.70mas1.01_APE0.4N50_smf0.2as0.12_s1": "M4",
    "chimera_out_64x32x32x32nc4nf2nas3b6.5mf0.72mas1.01_APE0.4N50_smf0.24as0.12_s1": "M5",
}

DATA = h5py.File("../tmp_data/chimera_data_full.hdf5")

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

DATA = h5py.File("../tmp_data/chimera_data_full.hdf5")


savename = "metadata/mesonF_meta.csv"
n = 0
if os.path.isfile(savename):
    meta = pd.read_csv(savename)

    existence = True

else:
    meta = []
    print("no", savename, "  creating a new one")
    existence = False


log_name = DATA.keys()

Nsource = "80"

for ens_tag in log_name:
    tmp = ens_tag.split("_")[2]
    ens = tmp.split("nc")[0] + "b" + tmp.split("b")[1]

    if existence and ens in meta.ENS.values:
        print(ens, "already done")
        continue

    print(ens)

    ens_group = DATA[ens_tag]

    Nt = int(ens.split("x")[0])
    Ns = int(ens.split("x")[0])
    beta = float(ens.split("b")[1].split("m")[0])
    f_bare_mass = ens_group["quarkmasses_fundamental"][0]
    epsilon_f = ens_group["Wuppertal_eps_fund"][0]

    data_tmp = []
    data_tmp.extend((ens, Nt, Ns, beta, f_bare_mass, epsilon_f, Nsource))

    for i in range(len(CHs)):
        ch = CHs[i]

        for Nsink in np.arange(0, 81, 10):
            tmp_bin = []
            for j in range(len(ch)):
                tmp_bin.append(
                    read_hdf.get_meson_corr(
                        DATA, ens_tag, "fund", Nsource, Nsink, ch[j]
                    )
                )

            corr = np.array(tmp_bin).mean(axis=0)

            c_boot = bootstrap.Correlator_resample(bootstrap.fold_correlators(corr))

            M_tmp = extract.Analysis_Mass_eff_cosh(
                c_boot,
                0,
                Nt,
                CHs_tag[i] + " [" + str(Nsource) + "," + str(Nsink) + "]",
            )

        plt.ylim(0.3, 0.5)
        plt.xlim(0, Nt / 2)
        extract.sperater.reset()

        plt.legend()
        plt.show()

        # N_source = smear.split('_')[1][1:]
        print(CHs_tag[i], end="  ")

        N_sink = input("N_sink = ")

        ti = int(input("ti = "))
        tf = int(input("tf = "))

        data_tmp.extend((N_sink, ti, tf))

    CSV_data.append(data_tmp)

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
            "N_source",
            "g5_N_sink",
            "g5_ti",
            "g5_tf",
            "gi_N_sink",
            "gi_ti",
            "gi_tf",
            "g0gi_N_sink",
            "g0gi_ti",
            "g0gi_tf",
            "g5gi_N_sink",
            "g5gi_ti",
            "g5gi_tf",
            "g0g5gi_N_sink",
            "g0g5gi_ti",
            "g0g5gi_tf",
            "id_N_sink",
            "id_ti",
            "id_tf",
        ]

        writer = csv.writer(csvfile)

        writer.writerow(fieldnames)
        writer.writerows(CSV_data)
