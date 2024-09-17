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
import plot_package


ens_tag = {
    "48x20x20x20b6.5mf0.71mas1.01": "chimera_out_48x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "64x20x20x20b6.5mf0.71mas1.01": "chimera_out_64x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "96x20x20x20b6.5mf0.71mas1.01": "chimera_out_96x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "64x20x20x20b6.5mf0.70mas1.01": "chimera_out_64x20x20x20nc4nf2nas3b6.5mf0.70mas1.01_APE0.4N50_smf0.2as0.12_s1",
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

CHs_tag = ["g5", "gi", "g0gi", "g5gi", "g0g5gi", "id"]
CHs_name = ["ps", "v", "t", "av", "at", "s"]


results = pd.read_csv("./CSVs/AS_meson_GEVP.csv")

for ens in list(results.ENS.values):
    print(ens)
    result_ens = results[results.ENS == ens]

    Nt = result_ens.Nt.values[0]

    for i in range(len(CHs)):
        ch = CHs[i]

        tmp_bin = []

        print(CHs_tag[i])

        for j in range(len(ch)):
            tmp_bin.append(
                read_hdf.get_meson_Cmat_single(
                    DATA, ens_tag[ens], "anti", 0, 80, 40, ch[j]
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
                M_tmp = extract.Analysis_Mass_eff_cosh(
                    LAM[:, :, n], 1, Nt / 2 + 2, f"E{n}"
                )
            else:
                E_string = fitting.print_non_zero(val, err)
                M_tmp = extract.Analysis_Mass_eff_cosh(
                    LAM[:, :, n], 1, Nt / 2 + 2, f"E{n} " + E_string
                )
                plot_package.plot_line(
                    val,
                    err,
                    result_ens.get(CHs_tag[i] + f"_E{n}_ti").values[0],
                    result_ens.get(CHs_tag[i] + f"_E{n}_tf").values[0],
                    plt.gca().lines[-1].get_color(),
                )

        extract.sperater.reset()
        plt.title(ens + " AS " + CHs_name[i])

        plt.ylim(0.3, 2)
        plt.legend(loc="upper right")
        plt.savefig("../plots/" + ens + "_AS_" + CHs_name[i] + ".pdf", transparent=True)
        plt.close()
        # plt.show()

results_mix = pd.read_csv("./CSVs/AS_meson_GEVP_mix.csv")

CHs_mix = [
    ["g1", "g0g1"],
    ["g2", "g0g2"],
    ["g3", "g0g3"],
]

for ens in list(results_mix.ENS.values):
    print(ens)
    result_ens = results_mix[results_mix.ENS == ens]

    Nt = result_ens.Nt.values[0]

    tmp_bin = []
    for i in range(len(CHs_mix)):
        ch = CHs_mix[i]
        tmp_bin.append(
            read_hdf.get_meson_Cmat_mix(
                DATA, ens_tag[ens], "anti", 0, 80, 40, ch[0], ch[1]
            )
        )

    Cmat = np.array(tmp_bin).mean(axis=0)

    t0_GEVP = result_ens.t0_GEVP.values[0]

    LAM_mix, VEC = extract.GEVP_fixT(Cmat, t0_GEVP, t0_GEVP + 1, Nt / 2 + 4)

    for n in range(Cmat.shape[-1]):
        val, err = (
            result_ens.get(f"vnt_E{n}").values[0],
            result_ens.get(f"vnt_E{n}_error").values[0],
        )

        if val == 0 and err == 0:
            M_tmp = extract.Analysis_Mass_eff_simple(
                LAM_mix[:, :, n], 1, Nt / 2 + 2, 1, f"E{n}"
            )

        else:
            E_string = fitting.print_non_zero(val, err)
            M_tmp = extract.Analysis_Mass_eff_simple(
                LAM_mix[:, :, n], 1, Nt / 2 + 2, 1, f"E{n} " + E_string
            )

            plot_package.plot_line(
                val,
                err,
                result_ens.get(f"vnt_E{n}_ti").values[0],
                result_ens.get(f"vnt_E{n}_tf").values[0],
                plt.gca().lines[-1].get_color(),
            )

    extract.sperater.reset()
    plt.title(ens + " AS v-t")

    plt.ylim(0.3, 2)
    plt.legend(loc="upper right")
    plt.savefig("../plots/" + ens + "_AS_vnt.pdf", transparent=True)
    plt.close()
