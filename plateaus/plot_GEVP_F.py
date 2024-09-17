import numpy as np
import h5py
import csv
import sys
import matplotlib.pyplot as plt
import os.path
import pandas as pd
import matplotlib.ticker as ticker

sys.path.insert(1, "./Lib")

import extract
import bootstrap
import read_hdf
import fitting
import plot_package


def LB_chimera(b):
    if b == "Lambda":
        return "Chimera_OC"
    elif b == "Sigma":
        return "Chimera_OV12"
    elif b == "SigmaS":
        return "Chimera_OV32"
    else:
        return "Chimera"


plt.figure(figsize=(7, 4.9))
plt.ylim(0.0, 2.0)
plt.tight_layout()
plt.title('')

ens_tag = {
    "48x20x20x20b6.5mf0.71mas1.01": "chimera_out_48x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "96x20x20x20b6.5mf0.71mas1.01": "chimera_out_96x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "64x20x20x20b6.5mf0.70mas1.01": "chimera_out_64x20x20x20nc4nf2nas3b6.5mf0.70mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "64x20x20x20b6.5mf0.71mas1.01": "chimera_out_64x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "64x32x32x32b6.5mf0.72mas1.01": "chimera_out_64x32x32x32nc4nf2nas3b6.5mf0.72mas1.01_APE0.4N50_smf0.24as0.12_s1",
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


CHs_name = ["PS", "V", "T", "AV", "AT", "S"]


result_mix = pd.read_csv("./CSVs/F_meson_GEVP_mix.csv")

CHs = [
    ["g1", "g0g1"],
    ["g2", "g0g2"],
    ["g3", "g0g3"],
]

CHs_tag = "VnT"

#ensembles = ['64x20x20x20b6.5mf0.71mas1.01','64x20x20x20b6.5mf0.70mas1.01','48x20x20x20b6.5mf0.71mas1.01','96x20x20x20b6.5mf0.71mas1.01','64x32x32x32b6.5mf0.72mas1.01']
ensembles = ['64x20x20x20b6.5mf0.71mas1.01']

for ens in ensembles:
    print(ens)
    result_ens = result_mix[result_mix.ENS == ens]

    Nt = result_ens.Nt.values[0]

    for i in range(len(CHs)):
        ch = CHs[i]

        tmp_bin = []
        plt.xlim(0.0, 28.5)
        print(ch)

        for j in range(len(ch)):
            tmp_bin.append(
                read_hdf.get_meson_Cmat_mix(
                    DATA, ens_tag[ens], "fund", 0, 80, 40, ch[0], ch[1]
                )
            )

    Cmat = np.array(tmp_bin).mean(axis=0)

    t0_GEVP = result_ens.t0_GEVP.values[0]

    LAM, VEC = extract.GEVP_fixT(Cmat, t0_GEVP, t0_GEVP + 1, Nt / 2 + 4)

    for n in range(Cmat.shape[-1]):
        val, err = (
            result_ens.get(CHs_tag + f"_E{n}").values[0],
            result_ens.get(CHs_tag + f"_E{n}_error").values[0],
        )

        if val == 0 and err == 0:
            M_tmp = extract.Analysis_Mass_eff_simple2(
                LAM[:, :, n], 1, Nt / 2 + 2, 1, f"E{n}", n, ' '
            )
        else:
            E_string = fitting.print_non_zero(val, err)
            M_tmp = extract.Analysis_Mass_eff_simple2(
                LAM[:, :, n], 1, Nt / 2 + 2, 1, f"E{n} " + E_string, n, E_string
            )

            plot_package.plot_line(
                val,
                err,
                result_ens.get(CHs_tag + f"_E{n}_ti").values[0],
                result_ens.get(CHs_tag + f"_E{n}_tf").values[0],
                plt.gca().lines[-1].get_color(),
            )

    extract.sperater.reset()
    
    
    # Set ticks every 5 units
    #plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(5))
    #plt.title(ens + " F " + CHs_tag)

    
    plt.legend(loc="upper right")
    plt.ylim(0.0, 2.0)
    plt.savefig("../plots/" + ens + "_F_" + CHs_tag + ".pdf", transparent=True)
    # plt.show()
    plt.close()




results = pd.read_csv("./CSVs/F_meson_GEVP.csv")

#ensembles = ['64x20x20x20b6.5mf0.71mas1.01','64x20x20x20b6.5mf0.70mas1.01','48x20x20x20b6.5mf0.71mas1.01','96x20x20x20b6.5mf0.71mas1.01','64x32x32x32b6.5mf0.72mas1.01']
ensembles = ['64x20x20x20b6.5mf0.71mas1.01']

CHs_tag = ["g5", "gi", "g0gi", "g5gi", "g0g5gi", "id"]

for ens in ensembles:
    print(ens)
    result_ens = results[results.ENS == ens]

    Nt = result_ens.Nt.values[0]

    for i in range(len(CHs)):
        ch = CHs[i]

        tmp_bin = []

        print(CHs_tag[i])
        
        plt.xlim(0.0, 24.5)
        
        for j in range(len(ch)):
            tmp_bin.append(
                read_hdf.get_meson_Cmat_single(
                    DATA, ens_tag[ens], "fund", 0, 80, 40, ch[j]
                )
            )

        Cmat = np.array(tmp_bin).mean(axis=0)

        t0_GEVP = result_ens.t0_GEVP.values[0]

        LAM, VEC = extract.GEVP_fixT(Cmat, t0_GEVP, t0_GEVP + 1, Nt / 2 + 4)

        for n in range(Cmat.shape[-1]):
            val, err = (
                result_ens.get(CHs_tag[i] + f"_E{n}").values[0],
                result_ens.get(CHs_tag[i] + f"_E{n}_error").values[0],
            )

            if val == 0 and err == 0:
                M_tmp = extract.Analysis_Mass_eff_cosh2(
                    LAM[:, :, n], 1, Nt / 2 + 2, f"E{n}", n, ' '
                )
            else:
                E_string = fitting.print_non_zero(val, err)
                M_tmp = extract.Analysis_Mass_eff_cosh2(
                    LAM[:, :, n], 1, Nt / 2 + 2, f"E{n} " + E_string, n, E_string
                )
                plot_package.plot_line(
                    val,
                    err,
                    result_ens.get(CHs_tag[i] + f"_E{n}_ti").values[0],
                    result_ens.get(CHs_tag[i] + f"_E{n}_tf").values[0],
                    plt.gca().lines[-1].get_color(),
                )

        extract.sperater.reset()
        #plt.title(ens + " F " + CHs_name[i])

        plt.ylim(0.0, 2)
        plt.legend(loc="upper right", fontsize=14)
        plt.savefig("../plots/" + ens + "_F_" + CHs_name[i] + ".pdf", transparent=True)
        plt.close()
        # plt.show()

