#!/usr/bin/env python
"""
Copyright (C) 2018-2024 Zidan Zhang <zhangzidan@gmail.com>

This file is distributed under free software licence:
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
from datetime import datetime
import time

import MDAnalysis as mda
import argparse
import logging
import os

import analysis_subh5
from analysis import *
from parameters import _args_analysis as _args

CHUNKS = 100


def subh5(h5, top, traj, anion, cation, begin, end, nt, cutoff, concentration, Dxyz):

    uta = mda.Universe(top, traj)

    h5tag_to, h5tag_bk, h5tag_if = 'total', 'bulk', 'interface'
    dtype = h5py.special_dtype(vlen=np.dtype('int, int'))
    dataset_to = get_or_create_dataset(h5, h5tag_to, dtype)
    dataset_bk = get_or_create_dataset(h5, h5tag_bk, dtype)
    dataset_if = get_or_create_dataset(h5, h5tag_if, dtype)

    frame_ids = [ts.frame for ts in uta.trajectory if (not begin or ts.time >= begin) and (not end or ts.time <= end)]

    do_analyse = functools.partial(analysis_subh5.analysis, top, traj, anion, cation, cutoff, concentration, Dxyz)

    pool = mp.Pool(processes=nt)

    data = pool.imap(do_analyse, frame_ids)

    # Aggregate data
    for _, tsHtt_to, tsHtt_bk, tsHtt_if in data:
        ds_append(dataset_to, tsHtt_to)
        ds_append(dataset_bk, tsHtt_bk)
        ds_append(dataset_if, tsHtt_if)


def aacf(h5filename, dataset_name, trestart, typeacf, nt):
    if typeacf == 'both':
        print('calculating C(t)')
        C_t = finalGetC_i(h5filename, dataset_name, 0, trestart, nt)
        acfCtPrint(dataset_name, C_t)
        print('calculating S(t)')
        S_t = finalGetC_c(h5filename, dataset_name, 0, trestart, nt)
        acfStPrint(dataset_name, S_t)
    elif typeacf == 'ct':
        print('calculating C(t)')
        C_t = finalGetC_i(h5filename, dataset_name, 0, trestart, nt)
        acfCtPrint(dataset_name, C_t)
    elif typeacf == 'st':
        print('calculating S(t)')
        S_t = finalGetC_c(h5filename, dataset_name, 0, trestart, nt)
        acfStPrint(dataset_name, S_t)


def main():
    args = _args().parse_args()
    logging.basicConfig(filename='utanalysis.log', level=logging.INFO, filemode='w')
    logging.info('STARTSTIME is: ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    if args.func == 'subh5':  # for homo polymer
        print('calculating subh5')
        if os.path.exists(args.h5):
            new_name = '{}_{}'.format(time.time(), args.h5)
            os.rename(args.h5, new_name)
            print(('renamed to {}_{}'.format(args.h5, new_name)))
        h5 = h5py.File(args.h5, 'w')
        subh5(h5, args.top, args.trj, args.anion, args.cation, args.begin, args.end, args.nt, args.cutoff, args.concentration, args.Dxyz)
        h5.close()
    elif args.func == 'aacf':  # for association time auto-correlation function
        print('calculating association auto-correlation function')
        h5 = h5py.File(args.h5, 'r')
        aacf(args.h5, args.h5tag, args.trestart, args.typeacf, args.nt)
        h5.close()
    logging.info('ENDSTIME is: ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


if __name__ == "__main__":
    main()
