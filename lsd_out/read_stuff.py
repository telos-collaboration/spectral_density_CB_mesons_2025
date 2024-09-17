import h5py

def print_hdf5_structure(group, path=''):
    for key in group.keys():
        item = group[key]
        new_path = f"{path}/{key}" if path else key
        if isinstance(item, h5py.Group):
            print_hdf5_structure(item, new_path)
        elif isinstance(item, h5py.Dataset):
            print(f"Dataset: {new_path} - Shape: {item.shape}, Dtype: {item.dtype}")

# Replace 'your_file.h5' with the path to your HDF5 file
file_path = 'chimera_baryon_data.hdf5'

try:
    with h5py.File(file_path, 'r') as file:
        print("HDF5 File Structure:")
        print_hdf5_structure(file)
except IOError as e:
    print(f"Error reading HDF5 file: {e}")
except Exception as e:
    print(f"An unexpected error occurred: {e}")

