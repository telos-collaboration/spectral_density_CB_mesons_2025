import h5py

def extract_dataset_column(file, dataset_path, column_index):
    try:
        with h5py.File(file, 'r') as hdf_file:
            dataset = hdf_file[dataset_path][()]
            column_data = dataset[:, column_index]
            return column_data
    except IOError as e:
        print(f"Error reading HDF5 file: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None


def extract_dataset(file, dataset_path, group1, group2):
    try:
        with h5py.File(file, 'r') as hdf_file:
            dataset = hdf_file[group1][group2][dataset_path][()]
            return dataset
    except IOError as e:
        print(f"Error reading HDF5 file: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None
