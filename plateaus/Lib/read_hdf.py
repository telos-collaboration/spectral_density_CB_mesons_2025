import h5py
import numpy as np
import sys

sys.path.insert(1, "./Lib")
import bootstrap
############# This is quite dependent on how Fabian structured the HDF5 file ########################

ensembles = ["M1", "M2", "M3", "M4", "M5"]
# Roots in HDF5 for each ensemble
roots = {
    "M1": "chimera_out_48x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "M2": "chimera_out_64x20x20x20b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "M3": "chimera_out_96x20x20x20nc4nf2nas3b6.5mf0.71mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "M4": "chimera_out_64x20x20x20nc4nf2nas3b6.5mf0.70mas1.01_APE0.4N50_smf0.2as0.12_s1",
    "M5": "chimera_out_64x32x32x32nc4nf2nas3b6.5mf0.72mas1.01_APE0.4N50_smf0.24as0.12_s1",
}

###### Functions to pic a specific correlator in a specific channel and ensemble from the HDF5 ######


def get_ensemble_root(ensemble_label):
    if ensemble_label in roots:
        return roots[ensemble_label]
    else:
        return None


def get_meson_corr(HDF5, ensemble, rep, Nsource, Nsink, channel):
    try:
        corr = HDF5[
            ensemble + "/" + f"source_N{Nsource}_sink_N{Nsink}/{rep} TRIPLET {channel}"
        ][0]

        return corr.T

    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None


"""
def get_meson_corr(file, ensemble , rep, Nsource, Nsink, channel):
    
    dataset_path = ensemble + "/" + f"source_N{Nsource}_sink_N{Nsink}/{rep} TRIPLET {channel}"
    
    try:
        with h5py.File(file, 'r') as hdf_file:
            dataset = hdf_file[dataset_path][0]
            
            return dataset.T
    except IOError as e:
        print(f"Error reading HDF5 file: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None
"""


def get_meson_Cmat_single(file, ensemble, rep, Nmin, Nmax, Nd, channel):
    source_Ns = [str(i) for i in np.arange(Nmin, Nmax + 1, Nd)]
    sink_Ns = [str(i) for i in np.arange(Nmin, Nmax + 1, Nd)]

    size = len(source_Ns)
    c_tmp = bootstrap.fold_correlators(
        get_meson_corr(file, ensemble, rep, Nmin, 0, channel)
    )

    mat = np.zeros(shape=(bootstrap.num_sample + 1, c_tmp.shape[1], size, size))

    # print(mat.shape)

    for i in range(size):
        for j in range(size):
            mat[:, :, i, j] = bootstrap.Correlator_resample(
                bootstrap.fold_correlators(
                    get_meson_corr(
                        file, ensemble, rep, source_Ns[i], sink_Ns[j], channel
                    )
                )
            )
    return mat


def get_meson_Cmat_mix(file, ensemble, rep, Nmin, Nmax, Nd, ch1, ch2):
    mixing_channel = [ch1, ch2]

    source_Ns = [str(i) for i in np.arange(Nmin, Nmax + 1, Nd)]
    sink_Ns = [str(i) for i in np.arange(Nmin, Nmax + 1, Nd)]
    size = len(source_Ns)

    c_tmp = bootstrap.Correlator_resample(
        get_meson_corr(file, ensemble, rep, Nmin, 0, ch1)
    )

    shape = list(c_tmp.shape)
    shape.append(size * len(mixing_channel))
    shape.append(size * len(mixing_channel))

    mat = np.zeros(shape=shape)

    for a in range(len(mixing_channel)):
        for b in range(len(mixing_channel)):
            for i in range(size):
                for j in range(size):
                    if a == b:
                        ch = mixing_channel[a]
                        mat[:, :, a * len(source_Ns) + i, b * len(sink_Ns) + j] = (
                            bootstrap.Correlator_resample(
                                bootstrap.fold_correlators(
                                    get_meson_corr(
                                        file,
                                        ensemble,
                                        rep,
                                        source_Ns[i],
                                        sink_Ns[j],
                                        ch,
                                    )
                                )
                            )
                        )

                    else:
                        ch = mixing_channel[a] + "_" + mixing_channel[b] + "_re"
                        mat[
                            :, :, a * len(source_Ns) + i, b * len(sink_Ns) + j
                        ] = -bootstrap.Correlator_resample(
                            bootstrap.fold_correlators_cross(
                                get_meson_corr(
                                    file, ensemble, rep, source_Ns[i], sink_Ns[j], ch
                                )
                            )
                        )

    return mat


########################### Example usage: ##########################################################
"""
file_path = "../tmp_data/correlators.hdf5"

# Replace these values with the actual values you want to use
Nsource_value = 80
Nsink_value = 40
rep ="fund"
channel = "g5"
ensemble = 'M3'   # M3 in this ensemble

result = bootstrap.Correlator_resample(get_meson_corr(file_path, ensemble, rep, Nsource_value, Nsink_value, channel))
"""
#####################################################################################################
