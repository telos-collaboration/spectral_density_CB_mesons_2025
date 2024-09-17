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
    ["g1", "g0g1"],
    ["g2", "g0g2"],
    ["g3", "g0g3"],
]


with open("metadata/MIX_GEVP_meta_AS.csv", newline="") as csvfile:
    meta = csv.DictReader(csvfile)
    n = 0
    for row in meta:
        print(
            30 * "~",
            row.get("ENS"),
            row.get("as_bare_mass"),
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
                row.get("as_bare_mass"),
                row.get("as_epsilon"),
                row.get("t0_GEVP"),
            )
        )

        t0_GEVP = int(row.get("t0_GEVP"))

        print(
            "---------------------------> GEVP with cross channels V & T <---------------------------"
        )

        tmp_bin = []
        for i in range(len(CHs)):
            ch = CHs[i]
            tmp_bin.append(
                read_hdf.get_meson_Cmat_mix(
                    DATA, ens_tag[row.get("ENS")], "anti", 0, 80, 40, ch[0], ch[1]
                )
            )

        Cmat = np.array(tmp_bin).mean(axis=0)

        LAM, VEC = extract.GEVP_fixT(Cmat, t0_GEVP, 1, float(row.get("Nt")))

        for en in range(Cmat.shape[-1]):
            print(
                "\n",
                15 * ">",
                f" E{en}",
                f"[Ti,Tf] = [%s,%s] " % (row.get(f"E{en}_ti"), row.get(f"E{en}_tf")),
                15 * "<",
                "\n",
            )

            if int(row.get(f"E{en}_ti")) == 0 and int(row.get(f"E{en}_tf")) == 0:
                m_tmp, m_tmp_err, chi2 = 0, 0, 0
                print("no plateau to fit...\n")

            else:
                m_tmp, m_tmp_err, chi2 = fitting.fit_exp_std(
                    LAM[:, :, en],
                    int(row.get(f"E{en}_ti")),
                    int(row.get(f"E{en}_tf")),
                    int(row.get("Nt")),
                )

            csv_row.extend(
                (
                    row.get(f"E{en}_ti"),
                    row.get(f"E{en}_tf"),
                    chi2,
                    m_tmp,
                    m_tmp_err,
                )
            )

        CSV_data_F.append(csv_row)

        n += 1


fieldnames = ["ENS", "Nt", "Ns", "beta", "as_bare_mass", "as_epsilon", "t0_GEVP"]


for n in range(Cmat.shape[-1]):
    fieldnames.extend(
        (
            f"vnt_E{n}_ti",
            f"vnt_E{n}_tf",
            f"vnt_E{n}_chisquare/dof",
            f"vnt_E{n}",
            f"vnt_E{n}_error",
        )
    )

with open("../CSVs/AS_meson_GEVP_mix.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(fieldnames)
    writer.writerows(CSV_data_F)
