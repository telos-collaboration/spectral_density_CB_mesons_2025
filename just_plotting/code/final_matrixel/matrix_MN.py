import argparse
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.legend_handler import HandlerTuple
import matplotlib.patches as mpatches
import h5py

parser = argparse.ArgumentParser()
parser.add_argument("--plot_styles", default="paperdraft.mplstyle")
parser.add_argument("--topology_h5", required=True)
args = parser.parse_args()

def read_specific_dataset(filename, dataset_path):
    with h5py.File(filename, 'r') as file:
        if dataset_path in file:
            return file[dataset_path][()]  # Read the full dataset
        else:
            raise KeyError(f"Dataset '{dataset_path}' not found in file '{filename}'")

# List of ensembles
ensembles = ['M1', 'M2', 'M3', 'M4', 'M5']

# Build dictionary
w0_values = {
    ens: read_specific_dataset(args.topology_h5, f"{ens}/w0_val")
    for ens in ensembles
}

files = [
    ('M1', ['input_fit/final_matrixel/M1_ground.txt']),
    ('M2', ['input_fit/final_matrixel/M2_ground.txt']),
    ('M3', ['input_fit/final_matrixel/M3_ground.txt']),
    ('M4', ['input_fit/final_matrixel/M4_ground.txt']),
    ('M5', ['input_fit/final_matrixel/M5_ground.txt']),
    ('M1', ['input_fit/final_matrixel/CB_M1_ground.txt']),
    ('M2', ['input_fit/final_matrixel/CB_M2_ground.txt']),
    ('M3', ['input_fit/final_matrixel/CB_M3_ground.txt']),
    ('M4', ['input_fit/final_matrixel/CB_M4_ground.txt']),
    ('M5', ['input_fit/final_matrixel/CB_M5_ground.txt'])
]

plt.style.use(args.plot_styles)
spacing = 0.
colors = [cm.tab10(i) for i in np.linspace(0, 0.5, 6)]
cb_colors = [cm.RdBu(i) for i in np.linspace(0, 0.3, 6)]
hatches = [r'//////', '++++']
offsets = {'M1': -0.05, 'M2': 0.0, 'M3': 0.05, 'M4': -0.25, 'M5': 0.50}
alpha_values = {'M1': 0.4, 'M2': 0.6, 'M3': 0.85, 'M4': 0.7, 'M5': 0.7}

x_labels = ['PS', 'V', 'AV', 'ps', 'v', 'av',  
            '$\Lambda^{+}_{\\rm CB}$', '$\Sigma^{+}_{\\rm CB}$', '$\Sigma^{*, +}_{\\rm CB}$', 
            '$\Lambda^{-}_{\\rm CB}$', '$\Sigma^{-}_{\\rm CB}$', '$\Sigma^{*, -}_{\\rm CB}$']
ciao = [0.3, 1.3, 2.3, 3.6, 4.6, 5.6, 7.1, 8.1, 9.1, 10.1, 11.1, 12.1]

def process_file(file_name, w0):
    values, errors = [], []
    with open(file_name, 'r') as file:
        for line in file:
            data = list(map(float, line.strip().split()))
            measurements, errors_data = data[::2], data[1::2]
            min_err_idx = errors_data.index(min(errors_data))
            stat = errors_data[min_err_idx]
            sys = max(abs(measurements[i] - measurements[j]) for i in range(len(measurements)) for j in range(i+1, len(measurements)))
            #print('ensemble: ', ensemble)
            #print('measurement: ', measurements[min_err_idx])
            #print('w0: ', w0)
            value = measurements[min_err_idx] * w0**2
            #print('value: ', value)
            total_err = math.sqrt(stat**2 + sys**2) * w0**2
            values.append(value)
            errors.append(total_err)
    return values, errors


def process_file2(file_name, w0):
    values, errors = [], []
    with open(file_name, 'r') as file:
        for line in file:
            data = list(map(float, line.strip().split()))
            measurements, errors_data = data[::2], data[1::2]
            min_err_idx = errors_data.index(min(errors_data))
            stat = errors_data[min_err_idx]
            sys = max(abs(measurements[i] - measurements[j]) for i in range(len(measurements)) for j in range(i+1, len(measurements)))
            #print('ensemble: ', ensemble)
            #print('measurement: ', measurements[min_err_idx])
            #print('w0: ', w0)
            value = measurements[min_err_idx] * w0**(3)
            #print('value: ', value)
            total_err = math.sqrt(stat**2 + sys**2) * w0**(3)
            values.append(value)
            errors.append(total_err)
    return values, errors

# ========== Plot 1: Fundamental + Antisymmetric ==========
fig1, ax1 = plt.subplots(figsize=(10, 6))
for ensemble, file_list in files:
    w0 = w0_values[ensemble]
    ensemble_offset = offsets[ensemble]

    for idx, file_name in enumerate(file_list):
        if 'CB' in file_name:
            continue  # Skip CB in first plot

        values, errors = process_file(file_name, w0)
        for i, (val, err) in enumerate(zip(values, errors)):
            if i == 3:
                err = 0.1*err
            x_pos = 0.3 + i + ensemble_offset
            if i >= 3:
                x_pos += 0.4
            if x_pos > 6.4:
                continue
            #print('Ensemble: ', ensemble)
            #print('val: ', val)   
            color = colors[i % 3]
            hatch = None
            hatch_color = None
            alpha = alpha_values.get(ensemble, 0.6)

            if ensemble == 'M4':
                hatch = hatches[0]
                hatch_color = colors[i % 3]
                color = 'none'
            elif ensemble == 'M5':
                hatch = hatches[1]
                hatch_color = colors[i % 3]
                color = 'none'

            rect = plt.Rectangle((x_pos - 0.2, val - err), 0.15 if ensemble in ['M4', 'M5'] else 0.4, 2 * err, color=color, alpha=alpha, hatch=hatch)
            if hatch_color:
                rect.set_edgecolor(hatch_color)
            ax1.add_patch(rect)

ax1.axvline(3.0, color='black', linestyle='-', alpha=0.7)
ax1.text(1.3, 0.75, 'Fundamental', fontsize=16, ha='center')
ax1.text(4.7, 0.75, 'Antisymmetric', fontsize=16, ha='center')
ax1.set_xticks(ciao[:6])
ax1.set_xticklabels(x_labels[:6], ha='right', fontsize=14, rotation=45)
ax1.set_ylabel('$\hat{c}_{M, 0}$', fontsize=16)
ax1.set_ylim(0.0, 0.8)
ax1.set_xlim(-0.5, 6.4)


import matplotlib.patches as mpatches

# Create custom legend handles with color patches for each group (same as your first code)
legend_handles = [
    mpatches.Patch(color='black', alpha=0.4),
    mpatches.Patch(color='black', alpha=0.6),
    mpatches.Patch(color='black', alpha=0.85),
    mpatches.Patch(color='black', alpha=0.7, fill=False),
    mpatches.Patch(color='black', alpha=0.7, fill=False),
]

hatches = [r'++++', '//////',]
for i, hatch in enumerate(hatches):
    legend_handles[-(i+1)].set_hatch(hatch)

legend_labels = ['M1', 'M2', 'M3', 'M4', 'M5']

# Add legend to the first plot
ax1.legend(
    handles=legend_handles,
    labels=legend_labels,
    handlelength=3.0,
    fontsize=14,
    loc='lower right',
    title_fontsize=14
)






plt.tight_layout()
plt.savefig('assets/plots/matrixel_fundamental_antisymmetric.pdf')
# plt.show()


# ========== Plot 2: Chimera Baryons ==========
fig2, ax2 = plt.subplots(figsize=(10, 6))
for ensemble, file_list in files:
    w0 = w0_values[ensemble]
    ensemble_offset = offsets[ensemble]

    for idx, file_name in enumerate(file_list):
        if 'CB' not in file_name:
            continue

        values, errors = process_file2(file_name, w0)
        for i, (val, err) in enumerate(zip(values, errors)):
            print('ensemble: ', ensemble)
            print('val: ', val)
            x_base = ciao[6 + i]
            x_pos = 0.3+ x_base + ensemble_offset

            color = cb_colors[i % 6]
            hatch = None
            hatch_color = None
            alpha = alpha_values.get(ensemble, 0.6)

            if ensemble == 'M4':
                hatch = hatches[1]
                hatch_color = cb_colors[i % 6]
                color = 'none'
            elif ensemble == 'M5':
                hatch = hatches[0]
                hatch_color = cb_colors[i % 6]
                color = 'none'

            rect = plt.Rectangle((x_pos - 0.2, val - err),
                                 0.15 if ensemble in ['M4', 'M5'] else 0.4,
                                 2 * err,
                                 color=color, alpha=alpha, hatch=hatch)
            if hatch_color:
                rect.set_edgecolor(hatch_color)
            ax2.add_patch(rect)

# Label
ax2.text(9.8, 0.09, 'Chimera baryons', fontsize=16, ha='center')

import matplotlib.patches as mpatches

# Create custom legend handles with color patches for each group (same as your first code)
legend_handles = [
    mpatches.Patch(color='black', alpha=0.4),
    mpatches.Patch(color='black', alpha=0.6),
    mpatches.Patch(color='black', alpha=0.85),
    mpatches.Patch(color='black', alpha=0.7, fill=False),
    mpatches.Patch(color='black', alpha=0.7, fill=False),
]

hatches = [r'++++', '//////']
for i, hatch in enumerate(hatches):
    legend_handles[-(i+1)].set_hatch(hatch)

legend_labels = ['M1', 'M2', 'M3', 'M4', 'M5']



# Add legend to the second plot
ax2.legend(
    handles=legend_handles,
    labels=legend_labels,
    handlelength=3.0,
    fontsize=14,
    loc='upper right',
    title_fontsize=14
)



ax2.set_xticks([x + 0.3 for x in ciao[6:]])
ax2.set_xticklabels(x_labels[6:], ha='right', fontsize=14, rotation=45)
ax2.set_ylabel('$\hat{K}_{B, 0}$', fontsize=16)
ax2.set_ylim(0.0, 0.1)
ax2.set_xlim(6.8, 13.0)
plt.tight_layout()
#plt.legend()
plt.savefig('assets/plots/matrixel_chimera_baryons.pdf')
