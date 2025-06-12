import argparse

import h5py


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_files", metavar="INPUT_FILE", nargs="+")
    parser.add_argument("--output_file", required=True)
    return parser.parse_args()


def main():
    args = get_args()
    with h5py.File(args.output_file, "w") as output_file:
        for input_filename in args.input_files:
            with h5py.File(input_filename, "r") as input_file:
                for key in input_file.keys():
                    if key in output_file:
                        raise ValueError(f"Key clash: {key}")
                    input_file.copy(key, output_file)


if __name__ == "__main__":
    main()
