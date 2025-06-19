import h5py


def extract_dataset_column(file, dataset_path, column_index):
    try:
        with h5py.File(file, "r") as hdf_file:
            dataset = hdf_file[dataset_path][()]
            # Check if dataset is 3D (e.g., shape (1, 48, 527)) or 2D (e.g., shape (48, 873))
            if len(dataset.shape) == 3 and dataset.shape[0] == 1:
                # Reshape (1, 48, 527) to (48, 527) for easier indexing
                dataset = dataset.reshape(dataset.shape[1], dataset.shape[2])
            # Extract the column from the 2D dataset
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
        with h5py.File(file, "r") as hdf_file:
            dataset = hdf_file[group1][group2][dataset_path][()]
            # Check if dataset is 3D (e.g., shape (1, 48, 527)) or 2D (e.g., shape (48, 873))
            if len(dataset.shape) == 3 and dataset.shape[0] == 1:
                # Reshape (1, 48, 527) to (48, 527)
                dataset = dataset.reshape(dataset.shape[1], dataset.shape[2])
            return dataset
    except IOError as e:
        print(f"Error reading HDF5 file: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None
