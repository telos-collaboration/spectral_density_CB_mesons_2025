import numpy as np
import json

def read_json(file_path):
    """Reads a JSON file and returns the parsed data."""
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

def bootstrap_errors(data, num_resamples=1000):
    """Calculates the bootstrap error from the data."""
    return np.mean(data), np.std(data)

def process_gevp_En_mass_samples(file_path, channel, rep, n):
    """Processes the 'gevp_f_at_E0_mass_samples' from the JSON file and calculates the average and bootstrap errors."""
    data = read_json(file_path)
    samples = np.array(data[f'gevp_{channel}_E{n}_mass_samples'])
    average, error = bootstrap_errors(samples)
    return average, error
