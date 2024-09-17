import numpy as np
import h5py
import csv
import sys
import matplotlib.pyplot as plt

sys.path.insert(1, "./Lib")

import extract
import bootstrap
import read_hdf


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


with open("metadata/mesonF_meta.csv", newline="") as csvfile:
    meta = csv.DictReader(csvfile)
    n = 0
    for row in meta:
        print(
            30 * "~",
            row.get("ENS"),
            row.get("f_bare_mass"),
            "Nsource=",
            row.get("N_source"),
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
                row.get("N_source"),
            )
        )

        for i in range(len(CHs)):
            ch = CHs[i]
            tmp_bin = []
            channel = CHs_tag[i]
            print(
                "-------> " + CHs_name[i],
                f": N_sink=%s [Ti,Tf] = [%d,%d] <-------"
                % (
                    row.get(channel + "_N_sink"),
                    int(row.get(channel + "_ti")),
                    int(row.get(channel + "_tf")),
                ),
            )
            for j in range(len(ch)):
                tmp_bin.append(
                    read_hdf.get_meson_corr(
                        DATA,
                        ens_tag[row.get("ENS")],
                        "fund",
                        row.get("N_source"),
                        row.get(channel + "_N_sink"),
                        ch[j],
                    )
                )

            corr = np.array(tmp_bin).mean(axis=0)
            bootstrap.num_confs = corr.shape[0]

            m_tmp, m_tmp_err, chi2 = extract.meson_mass(
                corr, int(row.get(channel + "_ti")), int(row.get(channel + "_tf"))
            )

            MASS[n, i] = m_tmp

            csv_row.extend(
                (
                    row.get(channel + "_N_sink"),
                    row.get(channel + "_ti"),
                    row.get(channel + "_tf"),
                    chi2,
                    m_tmp,
                    m_tmp_err,
                )
            )

        m_R = MASS[n, 0] / MASS[n, 1]
        m_R_err = bootstrap.bootstrap_error(m_R[0:-1], m_R[-1])

        csv_row.extend((m_R[-1], m_R_err))

        CSV_data_F.append(csv_row)

        n += 1

with open("../CSVs/F_meson.csv", "w", newline="") as csvfile:
    fieldnames = [
        "ENS",
        "Nt",
        "Ns",
        "beta",
        "f_bare_mass",
        "f_epsilon",
        "N_source",
        "PS_N_sink",
        "PS_ti",
        "PS_tf",
        "PS_chisquare/dof",
        "m_PS",
        "m_PS_error",
        "V_N_sink",
        "V_ti",
        "V_tf",
        "V_chisquare/dof",
        "m_V",
        "m_V_error",
        "T_N_sink",
        "T_ti",
        "V_tf",
        "T_chisquare/dof",
        "m_T",
        "m_T_error",
        "AV_N_sink",
        "AV_ti",
        "AV_tf",
        "AV_chisquare/dof",
        "m_AV",
        "m_AV_error",
        "AT_N_sink",
        "AT_ti",
        "AV_tf",
        "AT_chisquare/dof",
        "m_AT",
        "m_AT_error",
        "S_N_sink",
        "S_ti",
        "S_tf",
        "S_chisquare/dof",
        "m_S",
        "m_S_error",
        "m_PS/m_V",
        "m_PS/m_V_error",
    ]

    writer = csv.writer(csvfile)
    writer.writerow(fieldnames)
    writer.writerows(CSV_data_F)
