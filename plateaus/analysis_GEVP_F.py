import numpy as np
import h5py
import csv
import sys
import matplotlib.pyplot as plt

sys.path.insert(1, "./Lib")

import extract
import bootstrap
import read_hdf
import fitting


DATA = h5py.File("../input_correlators/chimera_data_full.hdf5")

CSV_data_F = []

MASS = np.zeros(shape=(5, 6, bootstrap.num_sample + 1))

ens_tag = {
    "48x20x20x20b6.5mf0.71mas1.01": "chimera_out_48x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "64x20x20x20b6.5mf0.71mas1.01": "chimera_out_64x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "96x20x20x20b6.5mf0.71mas1.01": "chimera_out_96x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "64x20x20x20b6.5mf0.70mas1.01": "chimera_out_64x20x20x20nc4nf2nas3b6.5mf0.70mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "64x32x32x32b6.5mf0.72mas1.01": "chimera_out_64x32x32x32nc4nf2nas3b6.5mf0.72mas1.01_APE0.4N50_smf0.24as0.12_s1",
}

CHs = [
    ["g5"],
    ["g1", "g2", "g3"],
    ["g0g1", "g0g2", "g0g3"],
    ["g5g1", "g5g2", "g5g3"],
    ["g0g5g1", "g0g5g2", "g0g5g3"],
    ["id"],
]

CHs_tag = ["g5", "gi", "g0gi", "g5gi", "g0g5gi", "id"]
CHs_name = ["PS", "V", "T", "AV", "AT", "S"]


with open("metadata/mesonF_GEVP_meta.csv", newline="") as csvfile:
    meta = csv.DictReader(csvfile)
    n = 0
    for row in meta:
        print(
            30 * "~",
            row.get("ENS"),
            row.get("f_bare_mass"),
            30 * "~",
            "\n",
        )

        csv_row = []

        csv_row.extend(
            (
                row.get("ENS"),
                row.get("Nt"),
                row.get("Ns"),
                row.get("beta"),
                row.get("f_bare_mass"),
                row.get("f_epsilon"),
                row.get("t0_GEVP"),
            )
        )

        for i in range(len(CHs)):
            ch = CHs[i]
            tmp_bin = []
            print(
                "--------------------------->",
                CHs_name[i],
                "<---------------------------",
            )
            for j in range(len(ch)):
                tmp_bin.append(
                    read_hdf.get_meson_Cmat_single(
                        DATA, ens_tag[row.get("ENS")], "fund", 0, 80, 40, ch[j]
                    )
                )

            Cmat = np.array(tmp_bin).mean(axis=0)

            t0_GEVP = int(row.get("t0_GEVP"))

            LAM, VEC = extract.GEVP_fixT(Cmat, t0_GEVP, 1, float(row.get("Nt")))

            for en in range(Cmat.shape[-1]):
                print(
                    "\n",
                    15 * ">",
                    f" E{en}",
                    f"[Ti,Tf] = [%s,%s] "
                    % (
                        row.get(CHs_tag[i] + f"_E{en}_ti"),
                        row.get(CHs_tag[i] + f"_E{en}_tf"),
                    ),
                    15 * "<",
                    "\n",
                )

                if (
                    int(row.get(CHs_tag[i] + f"_E{en}_ti")) == 0
                    and int(row.get(CHs_tag[i] + f"_E{en}_tf")) == 0
                ):
                    m_tmp, m_tmp_err, chi2 = 0, 0, 0
                    print("no plateau to fit...\n")

                else:
                    """
                    m_tmp, m_tmp_err, chi2 = fitting.X2_single_state_fit(
                        LAM[:, :, en],
                        int(row.get(CHs_tag[i] + f"_E{en}_ti")),
                        int(row.get(CHs_tag[i] + f"_E{en}_tf")),
                    )
                    """

                    m_tmp, m_tmp_err, chi2 = fitting.fit_cosh_std(
                        LAM[:, :, en],
                        int(row.get(CHs_tag[i] + f"_E{en}_ti")),
                        int(row.get(CHs_tag[i] + f"_E{en}_tf")),
                        int(row.get("Nt")),
                    )

                    # MASS[n, i ] = m_tmp

                csv_row.extend(
                    (
                        row.get(CHs_tag[i] + f"_E{en}_ti"),
                        row.get(CHs_tag[i] + f"_E{en}_tf"),
                        chi2,
                        m_tmp,
                        m_tmp_err,
                    )
                )

        CSV_data_F.append(csv_row)

        n += 1

np.save("tmp_data/GEVP_PS_F.npy", MASS)

fieldnames = ["ENS", "Nt", "Ns", "beta", "f_bare_mass", "f_epsilon", "t0_GEVP"]

for ch in CHs_tag:
    for n in range(Cmat.shape[-1]):
        fieldnames.extend(
            (
                f"{ch}_E{n}_ti",
                f"{ch}_E{n}_tf",
                f"{ch}_E{n}_chisquare/dof",
                f"{ch}_E{n}",
                f"{ch}_E{n}_error",
            )
        )


with open("../CSVs/F_meson_GEVP.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(fieldnames)
    writer.writerows(CSV_data_F)
