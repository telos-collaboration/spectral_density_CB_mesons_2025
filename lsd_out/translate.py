import numpy as np


def save_matrix_to_file(matrix, output_file):
    N, R = matrix.shape[2], matrix.shape[1]

    with open(output_file, 'w') as output_file:
        output_file.write(f'{N} {R} {R} {2} {3}\n')
        for col in range(N):
            for row in range(R):
                output_file.write(f'{row % R} {matrix[0, row, col]}\n')
