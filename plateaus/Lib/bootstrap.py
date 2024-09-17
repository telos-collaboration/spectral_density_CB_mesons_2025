import numpy as np
import random

### hard coded bootstrap parameters

bin_size = 2
num_sample = 1000
num_confs = 0
seed = 8437

######


def random_index(num_confs):
    random.seed(seed)

    random_index = np.zeros(shape=(num_sample, int(num_confs / bin_size)), dtype=int)
    for i in range(num_sample):
        for j in range(int(num_confs / bin_size)):
            random_index[i, j] = random.randint(0, int(num_confs / bin_size) - 1)

    return random_index


def bootstrap_main(C):
    binS = []
    num_confs = C.shape[0]

    num_bin = int(num_confs / bin_size)

    index_list = random_index(num_confs)

    for i in range(num_bin):
        b = C[bin_size * i : bin_size * (i + 1), :]
        binS.append(np.mean(b, axis=0))

    binS = np.array(binS)

    bootstrap_set = []
    for i in range(num_sample):
        sample = []
        for j in range(num_bin):
            sample.append(binS[index_list[i, j]])

        bootstrap_set.append(np.mean(sample, axis=0))

    return np.array(bootstrap_set)


def bootstrap_error(A, centval):
    B = np.copy(A)
    B.sort()
    index_err = int(len(B) * 0.18)

    MEAN = np.mean(B)
    ERROR = max(abs(B[-index_err] - centval), abs(centval - B[index_err - 1]))

    return ERROR


def fold_correlators(C):
    return (C + np.roll(np.flip(C, axis=1), 1, axis=1)) / 2


def fold_correlators_cross(C):
    C_fold = (C - np.roll(np.flip(C, axis=1), 1, axis=1)) / 2

    C_fold[:, 0] = C[:, 0]

    return C_fold


def gauB_DIS(cen, err, size):
    atmp = np.random.normal(cen, scale=err, size=size)
    return np.append(atmp, cen)


def Correlator_resample(C_ALL):
    C_ALL_resample = np.zeros(shape=(num_sample + 1, np.shape(C_ALL)[1]))
    C_ALL_resample[0:num_sample] = bootstrap_main(C_ALL)
    C_ALL_resample[num_sample] = np.mean(C_ALL, axis=0)

    return C_ALL_resample
