import argparse
import json
import numpy as np
import pandas as pd
import os
import re
import h5py
import math

# Constants
PI2 = math.pi**2

parser = argparse.ArgumentParser()
parser.add_argument("--topology_h5", required=True)
args = parser.parse_args()

# Read the renormalise.csv file
df = pd.read_csv("./metadata/renormalise.csv")

# Compute Z values dynamically using CSV and HDF5 plaquette data
computed_z = {}

with h5py.File(args.topology_h5, "r") as f:
    for _, row in df.iterrows():
        ens = row["Ens"]
        beta = row["beta"]
        C_fund = row["C_fund"]
        C_as = row["C_as"]
        Delta_Sigma1 = row["Delta_Sigma1"]
        Delta_R1 = row["Delta_R1"]
        Delta_R2 = row["Delta_R2"]
        Delta_Lambda = row["Delta_Lambda"]
        Delta_Sigma = row["Delta_Sigma"]

        # Load and average plaquette
        plaq_path = f"{ens}/plaquette"
        if plaq_path not in f:
            raise KeyError(f"Missing {plaq_path} in {args.topology_h5}")
        plaq_values = f[plaq_path][()]
        plaq_avg = np.mean(plaq_values)

        # Compute Zs
        factor_fund = 8 * C_fund / (16 * PI2 * beta * plaq_avg)
        factor_as = 8 * C_as / (16 * PI2 * beta * plaq_avg)

        computed_z[ens] = {
            "Z_PS_fund": 1 + factor_fund * (Delta_Sigma1 + Delta_R1),
            "Z_V_fund": 1 + factor_fund * (Delta_Sigma1 + Delta_R2),
            "Z_A_fund": 1 + factor_fund * (Delta_Sigma1 + Delta_R1),
            "Z_PS_as": 1 + factor_as * (Delta_Sigma1 + Delta_R1),
            "Z_V_as": 1 + factor_as * (Delta_Sigma1 + Delta_R2),
            "Z_A_as": 1 + factor_as * (Delta_Sigma1 + Delta_R1),
            "Z_Lambda": 1
            + (factor_fund / C_fund)
            * ((C_fund + 1 / 2 * C_as) * Delta_Sigma1 + Delta_Lambda),
            "Z_Sigma": 1
            + (factor_fund / C_fund)
            * ((C_fund + 1 / 2 * C_as) * Delta_Sigma1 + Delta_Sigma),
        }

print(computed_z)


# Ensemble directory to renormalisation key mapping
ensemble_map = {
    "./JSONs/Sp4b6.5nF2nAS3mF-0.71mAS-1.01T48L20": "M1",
    "./JSONs/Sp4b6.5nF2nAS3mF-0.71mAS-1.01T64L20": "M2",
    "./JSONs/Sp4b6.5nF2nAS3mF-0.71mAS-1.01T96L20": "M3",
    "./JSONs/Sp4b6.5nF2nAS3mF-0.7mAS-1.01T64L20": "M4",
    "./JSONs/Sp4b6.5nF2nAS3mF-0.72mAS-1.01T64L32": "M5",
}

# Load renormalization factors
renorm_df = pd.read_csv("./metadata/renormalise.csv")

# Results will be stored here
results = []

for ensemble_dir, ens_key in ensemble_map.items():
    print(f"\n==============================")
    print(f"Ensemble: {ensemble_dir} (Ens = {ens_key})")

    try:
        # --- Load fundamental JSONs ---
        with open(
            os.path.join(ensemble_dir, "meson_extraction_f_ps_samples.json")
        ) as f:
            f_ps_data = json.load(f)
        with open(
            os.path.join(ensemble_dir, "meson_extraction_f_av_samples.json")
        ) as f:
            f_av_data = json.load(f)
        with open(os.path.join(ensemble_dir, "meson_extraction_f_v_samples.json")) as f:
            f_v_data = json.load(f)

        # --- Load antisymmetric JSONs ---
        with open(
            os.path.join(ensemble_dir, "meson_extraction_as_ps_samples.json")
        ) as f:
            as_ps_data = json.load(f)
        with open(
            os.path.join(ensemble_dir, "meson_extraction_as_av_samples.json")
        ) as f:
            as_av_data = json.load(f)
        with open(
            os.path.join(ensemble_dir, "meson_extraction_as_v_samples.json")
        ) as f:
            as_v_data = json.load(f)
    except FileNotFoundError as e:
        print(f"Missing file in {ensemble_dir}: {e}")
        continue

    try:
        renorm_row = computed_z[ens]
    except IndexError:
        print(f"No renormalization data for Ens = {ens_key}")
        continue

    def compute_observables(m_ps, m_av, m_v, me_ps, me_av, me_v, Z_PS, Z_A, Z_V):
        f_ps = me_ps / (np.sqrt(2) * m_ps) * Z_PS
        f_av = me_av / (np.sqrt(2) * m_av) * Z_A
        f_v = me_v / (np.sqrt(2) * m_v) * Z_V

        s_0 = 4 * np.pi * ((f_v**2 / m_v**2) - (f_av**2 / m_av**2))
        s_1 = 1 - ((f_av**2 + f_ps**2) / f_v**2)
        s_2 = 1 - ((m_av**2 * f_av**2) / (m_v**2 * f_v**2))

        return {
            "s_0": (np.mean(s_0), np.std(s_0, ddof=1)),
            "s_1": (np.mean(s_1), np.std(s_1, ddof=1)),
            "s_2": (np.mean(s_2), np.std(s_2, ddof=1)),
        }

    # Fundamental
    obs_f = compute_observables(
        np.array(f_ps_data["f_ps_mass_samples"]),
        np.array(f_av_data["f_av_mass_samples"]),
        np.array(f_v_data["f_v_mass_samples"]),
        np.array(f_ps_data["f_ps_matrix_element"]),
        np.array(f_av_data["f_av_matrix_element"]),
        np.array(f_v_data["f_v_matrix_element"]),
        renorm_row["Z_PS_fund"],
        renorm_row["Z_A_fund"],
        renorm_row["Z_V_fund"],
    )
    results.append(
        {
            "representation": "f",
            "Ens": ens_key,
            "s_0": obs_f["s_0"][0],
            "s_0_err": obs_f["s_0"][1],
            "s_1": obs_f["s_1"][0],
            "s_1_err": obs_f["s_1"][1],
            "s_2": obs_f["s_2"][0],
            "s_2_err": obs_f["s_2"][1],
        }
    )

    # Antisymmetric
    obs_as = compute_observables(
        np.array(as_ps_data["as_ps_mass_samples"]),
        np.array(as_av_data["as_av_mass_samples"]),
        np.array(as_v_data["as_v_mass_samples"]),
        np.array(as_ps_data["as_ps_matrix_element"]),
        np.array(as_av_data["as_av_matrix_element"]),
        np.array(as_v_data["as_v_matrix_element"]),
        renorm_row["Z_PS_as"],
        renorm_row["Z_A_as"],
        renorm_row["Z_V_as"],
    )
    results.append(
        {
            "representation": "as",
            "Ens": ens_key,
            "s_0": obs_as["s_0"][0],
            "s_0_err": obs_as["s_0"][1],
            "s_1": obs_as["s_1"][0],
            "s_1_err": obs_as["s_1"][1],
            "s_2": obs_as["s_2"][0],
            "s_2_err": obs_as["s_2"][1],
        }
    )

# Convert results to DataFrame and save to CSV
results_df = pd.DataFrame(results)
results_df = results_df.sort_values(by=["representation", "Ens"])
results_df.to_csv("./CSVs/s_parameters_summary.csv", index=False)

# Generate LaTeX-style table
print("\n=== LaTeX Table ===")
for rep in ["f", "as"]:
    for ens in ["M1", "M2", "M3", "M4", "M5"]:
        row = results_df[
            (results_df["representation"] == rep) & (results_df["Ens"] == ens)
        ]
        if row.empty:
            continue
        r = row.iloc[0]

        def fmt(val, err):
            return f"{val:.3f}({int(round(err * 1000)):d})"

        print(
            f"$({{{rep}}})$-type     & {ens} & {fmt(r['s_0'], r['s_0_err'])} & {fmt(r['s_1'], r['s_1_err'])} & {fmt(r['s_2'], r['s_2_err'])} \\\\"
        )
    if rep == "f":
        print("        \\hline")


# Write LaTeX table to file
with open("./assets/tables/s_parameters_table.tex", "w") as texfile:
    texfile.write("\\begin{tabular}{ |c|c|c|c|c| }\n")
    texfile.write("    \\hline\\hline\n")
    texfile.write("    Meson  & Ensemble & $s_0$ & $s_1$ & $s_2$ \\\\\n")
    texfile.write("    \\hline\n")

    for rep in ["f", "as"]:
        for ens in ["M1", "M2", "M3", "M4", "M5"]:
            row = results_df[
                (results_df["representation"] == rep) & (results_df["Ens"] == ens)
            ]
            if row.empty:
                continue
            r = row.iloc[0]

            def fmt(val, err):
                return f"{val:.3f}({int(round(err * 1000)):d})"

            texfile.write(
                f"$({{\\rm {rep}}})$-type     & {ens} & {fmt(r['s_0'], r['s_0_err'])} & {fmt(r['s_1'], r['s_1_err'])} & {fmt(r['s_2'], r['s_2_err'])} \\\\\n"
            )
        if rep == "f":
            texfile.write("    \\hline\n")

    texfile.write("    \\hline  \\hline\n")
    texfile.write("\\end{tabular}\n")
