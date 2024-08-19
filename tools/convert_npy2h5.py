#! /usr/bin/env python

import argparse
import h5py
import numpy as np


args = argparse.ArgumentParser()
args.add_argument('--npy', required=True)
args.add_argument('--dataset', default=None)
args.add_argument('--h5', required=True)

argv = args.parse_args()

data = np.load(argv.npy)

h5 = h5py.File(argv.h5, 'a')
if argv.dataset not in h5:
    print('Dataset {} not found in {}. Create'.format(argv.dataset, argv.h5))
    dtype = h5py.special_dtype(vlen=np.dtype('int, int'))
    dataset = h5.create_dataset(argv.dataset, (len(data), ), dtype=dtype)
else:
    dataset = h5[argv.dataset]

for i in range(data.shape[0]):
    dataset[i] = data[i]

h5.close()
print('Moved {} to {}'.format(argv.npy, argv.h5))
