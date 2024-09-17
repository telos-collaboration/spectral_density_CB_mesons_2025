import numpy as np
import h5py
import csv
import sys
import matplotlib.pyplot as plt
import os.path
import pandas as pd

sys.path.insert(1, "./Lib")

import read_hdf
import plot_package


ens_tag = {
    "48x20x20x20b6.5mf0.71mas1.01": "M1",
    "64x20x20x20b6.5mf0.71mas1.01": "M2",
    "96x20x20x20b6.5mf0.71mas1.01": "M3",
    "64x20x20x20b6.5mf0.70mas1.01": "M4",
    "64x32x32x32b6.5mf0.72mas1.01": "M5",
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
CHs_name = {"g5": "PS", "gi": "V", "g0gi": "T", "g5gi": "AV", "g0g5gi": "AT", "id": "S"}


results = pd.read_csv("CSVs/F_meson_GEVP.csv")

for ch in CHs_tag:
    plt.errorbar(
        np.arange(0, 5),
        results.get(f"{ch}_E0").values,
        results.get(f"{ch}_E0_error").values,
        marker="o",
        linestyle="",
        label=CHs_name[ch],
        alpha=0.6,
    )

ticks = []
for i in range(5):
    ticks.append(ens_tag[results.ENS.values[i]])

plt.xticks(np.arange(0, 5), ticks)
plt.legend(loc="upper right")
plt.show()
