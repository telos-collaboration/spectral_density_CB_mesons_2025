import h5py

def get_root_paths(h5_file):
    root_paths = set()

    def visit(name):
        root_paths.add(name.split('/')[0])

    with h5py.File(h5_file, 'r') as file:
        file.visit(visit)

    return root_paths

hdf5_file = 'chimera_baryon_data.hdf5'
roots = get_root_paths(hdf5_file)

print("Unique root paths in the HDF5 file:")
for root in roots:
    print(root)

