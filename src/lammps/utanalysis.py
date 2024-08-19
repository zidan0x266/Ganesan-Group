#!/usr/bin/env python
"""
Copyright (C) 2018-2021 Zidan Zhang <zhangzidan@gmail.com>
Jakub Krajniak <jkrajniak@gmail.com>

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

import analysis_copoly
import analysis_homopoly
import analysis_hopping
from analysis import *
from parameters import _args_analysis as _args

CHUNKS = 100


def copoly(h5, h5tag, top, traj, anion, cation, beadiny, begin, end, nt, cutoff, Icutoff, counterN):
    """
    This function outputs the analysis for everything of interest for
    copolymer (random or block), including association relationship,
    time auto correlation function and these properties near the
    interface.
    """
    
    uta = mda.Universe(top, traj, format="LAMMPS")
    # Store the associations for 1. anion and cation (AC); 2. anion and chain (CH)
    freqAC = []
    freqCH = []
    freqBAC = []
    freqBCH = []  # at bulk
    freqIAC = []
    freqICH = []  # at interface
    anionIf = []  # Anions at interface

    in_h5tag = 'in_{}'.format(h5tag)
    it_h5tag = 'it_{}'.format(h5tag)

    dtype = h5py.special_dtype(vlen=np.dtype('int, int'))
    ds_htt = get_or_create_dataset(h5, h5tag, dtype)
    ds_in_htt = get_or_create_dataset(h5, in_h5tag, dtype)
    ds_it_htt = get_or_create_dataset(h5, it_h5tag, dtype)

    frame_ids = [ts.frame for ts in uta.trajectory if (not begin or ts.time >= begin) and (not end or ts.time <= end)]

    do_analyse = functools.partial(analysis_copoly.analysis, top, traj, anion, cation, beadiny, cutoff, Icutoff, counterN)

    pool = mp.Pool(processes=nt)

    data = pool.imap(do_analyse, frame_ids)

    for _, tsAC, tsCH, tsBAC, tsBCH, tsIAC, tsICH, tsIAE, tsHtt, tsInHtt, tsItHtt in data:
        freqAC.append(tsAC)
        freqCH.append(tsCH)
        freqBAC.append(tsBAC)
        freqBCH.append(tsBCH)
        freqIAC.append(tsIAC)
        freqICH.append(tsICH)
        anionIf.append(tsIAE)

        ds_append(ds_htt, tsHtt)
        ds_append(ds_in_htt, tsInHtt)
        ds_append(ds_it_htt, tsItHtt)

    return freqAC, freqCH, freqBAC, freqBCH, freqIAC, freqICH, anionIf


def homopoly(h5, h5tag, top, traj, anion, cation, begin, end, nt, cutoff, counterN):
    """
    This function outputs the analysis for everything of interest for
    homopolymer, including association relationship, time auto 
    correlation function.
    """

    uta = mda.Universe(top, traj, format="LAMMPS")
    # Store the associations for 1. anion and cation (AC); 2. anion and chain (CH)
    freqAC = []
    freqCH = []

    dtype = h5py.special_dtype(vlen=np.dtype('int, int'))
    dataset = get_or_create_dataset(h5, h5tag, dtype)

    frame_ids = [ts.frame for ts in uta.trajectory if (not begin or ts.time >= begin) and (not end or ts.time <= end)]

    do_analyse = functools.partial(analysis_homopoly.analysis, top, traj, anion, cation, cutoff, counterN)

    pool = mp.Pool(processes=nt)

    data = pool.imap(do_analyse, frame_ids)

    # Aggregate data
    for _, tsAC, tsCH, tsHtt in data:
        freqAC.append(tsAC)
        freqCH.append(tsCH)
        ds_append(dataset, tsHtt)

    return freqAC, freqCH


def hopping(h5, dataset_name, counterN):
    hopping, dissociation, association = analysis_hopping.analysis(h5, dataset_name, counterN)
    hopesPrint('hopping', hopping[0], hopping[1])
    hopesPrint('dissociation', dissociation[0], dissociation[1])
    hopesPrint('association', association[0], association[1])


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

    if args.func == 'copoly':  # for random or block copolymer
        print('calculating copoly')
        if os.path.exists(args.h5):
            new_name = '{}_{}'.format(time.time(), args.h5)
            os.rename(args.h5, new_name)
            print(('renamed to {}_{}'.format(args.h5, new_name)))
        h5 = h5py.File(args.h5, 'w')
        freqAC, freqCH, freqBAC, freqBCH, freqIAC, freqICH, anionIf = copoly(h5,
                                                           args.h5tag,
                                                           args.top,
                                                           args.trj,
                                                           args.anion,
                                                           args.cation,
                                                           args.beadiny,
                                                           args.begin,
                                                           args.end,
                                                           args.nt,
                                                           args.cutoff,
                                                           args.Icutoff,
                                                           args.counterN)
        acNum, acPro = assoFreq(freqAC)
        chNum, chPro = assoFreq(freqCH)
        acBNum, acBPro = assoFreq(freqBAC)
        chBNum, chBPro = assoFreq(freqBCH)
        acINum, acIPro = assoFreq(freqIAC)
        chINum, chIPro = assoFreq(freqICH)
        anaPrint('assoPair', acNum, acPro)
        anaPrint('assoChain', chNum, chPro)
        anaPrint('assoBPair', acBNum, acBPro)
        anaPrint('assoBChain', chBNum, chBPro)
        anaPrint('assoIPair', acINum, acIPro)
        anaPrint('assoIChain', chINum, chIPro)
        interFPrint('anionIf', np.average(anionIf), np.std(anionIf))
        h5.close()
    elif args.func == 'homopoly':  # for homo polymer
        print('calculating homopoly')
        if os.path.exists(args.h5):
            new_name = '{}_{}'.format(time.time(), args.h5)
            os.rename(args.h5, new_name)
            print(('renamed to {}_{}'.format(args.h5, new_name)))
        h5 = h5py.File(args.h5, 'w')
        freqAC, freqCH = homopoly(h5, args.h5tag, args.top, args.trj, args.anion, args.cation, 
                                    args.begin, args.end, args.nt, args.cutoff, args.counterN)
        # rawSave(h5, 'htt', ISAACS)
        acNum, acPro = assoFreq(freqAC)
        chNum, chPro = assoFreq(freqCH)
        anaPrint('assoPair', acNum, acPro)
        anaPrint('assoChain', chNum, chPro)
        h5.close()
    elif args.func == 'hopping':  # for hopping event analysis
        print('calculating hopping events')
        h5 = h5py.File(args.h5, 'r')
        hopping(h5, args.h5tag, args.counterN)
        h5.close()
    elif args.func == 'aacf':  # for association time auto-correlation function
        print('calculating association auto-correlation function')
        h5 = h5py.File(args.h5, 'r')
        aacf(args.h5, args.h5tag, args.trestart, args.typeacf, args.nt)
        h5.close()
    logging.info('ENDSTIME is: ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S'))


if __name__ == "__main__":
    main()
