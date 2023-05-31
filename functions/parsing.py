import h5py
import numpy as np
import bz2
import _pickle as cPickle
import pickle
import os
from collections import defaultdict
from tqdm import tqdm
from glob import glob


def mat2py(mat_path):
    '''This function converts matlab structs to python dictionary (files names are keys and values are filed content).'''
    with h5py.File(mat_path) as f:
        keys = list(f.keys())
        mat_struct = f[keys[1]]
        field_dict = defaultdict(list)
        for field_name in tqdm(mat_struct):
            for ind in range(len(mat_struct[field_name])):
                field_dict[field_name].append(np.array(f[mat_struct[field_name][ind, 0]]))
    return field_dict


def compressed_pickle(title, data):
    '''Pickle a file and then compress it into a file with extension.'''
    with bz2.BZ2File(title + '.pbz2', 'w') as f:
        cPickle.dump(data, f)


def decompress_pickle(file):
    '''Load any compressed pickle file.'''
    data = bz2.BZ2File(file, 'rb')
    data = cPickle.load(data)
    return data


def load_dict(path):
    fin_dict = dict()
    files = glob(os.path.join(path, '*.pkl'))
    for file in files:
        with open(file, 'rb') as f:
            fin_dict.update(pickle.load(f))
    return fin_dict


if __name__ == '__main__':
    x = 1