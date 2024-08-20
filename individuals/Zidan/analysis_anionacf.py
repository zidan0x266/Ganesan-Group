"""
Copyright (C) 2018-2020 Zidan Zhang <zhangzidan@gmail.com>
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
import multiprocessing as mp
import functools

import MDAnalysis as mda
import numpy as np
import h5py


def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def get_or_create_dataset(h5, dataset_name, dtype):
    if dataset_name not in h5:
        return h5.create_dataset(dataset_name, (0, ), maxshape=(None, ), dtype=dtype, chunks=(100, ))
    return h5[dataset_name]


def ds_append(ds, val):
    ds.resize(ds.shape[0]+1, axis=0)
    ds[-1] = val


def anionPrint(filename, acf):
    time = np.linspace(0, len(acf) - 1, len(acf))
    with open('{}.xvg'.format(filename), 'w') as sqtout:
        print('# time, auto-correlation function', file=sqtout)
        for i in range(len(time)):
            print('{:10.3f} {:10.5f}'.format(time[i], acf[i]), file=sqtout)


def getC_i(AGS, t0, t):
    """
    This function give the intermitent Association Lifetime
    C_i = <h(t0)h(t)>/<h(t0)> between t0 and t
    """
    if AGS[t0].shape[0] == 0:
        return 0.0

    Hit = np.intersect1d(AGS[t0], AGS[t])
    C_i = len(Hit) / float(len(AGS[t0]))

    return C_i


def intervC_i(h5filename, dataset_name, t0, tf, reservoir, dt):
    """
    This function gets all the data for the h(t0)h(t0+dt), where
    t0 = 1,2,3,...,tf. This function give us a point of the final plot
    C(t) vs t,
    """
    AGS, h5 = rawLoad(h5filename, dataset_name)
    a = 0
    count = 0
    for t0 in reservoir:
        if t0 + dt <= tf:
            a += getC_i(AGS, t0, t0 + dt)
            count += 1
    h5.close()
    return a / count


def finalGetC_i(h5filename, dataset_name, t0, trestart, nt):
    """
    This function gets the final data of the C_i graph.
    """
    AGS, h5 = rawLoad(h5filename, dataset_name)
    tf = AGS.shape[0] - 1
    num_AGS = AGS.shape[0]
    h5.close()

    pool = mp.Pool(processes=nt)
    reservoir = [x * trestart for x in range(tf) if x * trestart < tf]
    func = functools.partial(intervC_i, h5filename, dataset_name, t0, tf, reservoir)
    return_data = pool.map(func, list(range(num_AGS)))
    return return_data


def getC_c(AGS, t0, t):
    if t0 == t:
        return 1.0

    dt = 1
    ht0 = AGS[t0]
    htt_old = AGS[t0]

    while t0 + dt <= t:
        htt = AGS[t0 + dt]
        Htt = np.intersect1d(htt_old, htt)
        t0 += dt
        htt_old = Htt
    C_c = len(Htt) / float(len(ht0))
    if AGS[t0].shape[0] == 0:
        return 0
    else:
        return C_c


def intervC_c(h5filename, dataset_name, t0, tf, reservoir, dt):
    AGS, h5 = rawLoad(h5filename, dataset_name)
    a = 0
    count = 0
    for t0 in reservoir:
        if t0 + dt <= tf:
            if t0 == t0 + dt:
                break
            a += getC_c(AGS, t0, t0 + dt)
            count += 1

    h5.close()

    if count == 0:
        return 1.0

    return a / count


def finalGetC_c(h5filename, dataset_name, t0, trestart, nt):
    AGS, h5 = rawLoad(h5filename, dataset_name)
    tf = AGS.shape[0] - 1
    num_AGS = AGS.shape[0]
    h5.close()

    pool = mp.Pool(processes=nt)
    reservoir = [x * trestart for x in range(tf) if x * trestart < tf]
    func = functools.partial(intervC_c, h5filename, dataset_name, t0, tf, reservoir)
    return_data = pool.map(func, list(range(num_AGS)))
    return return_data


def getgroups(frame, ASSOac, ASSOal):
    at1, at2 = [], []
    anion = list(set([x[0] for x in ASSOal[frame]]))  # find the anions in co-coordination
    for pair in ASSOac[frame]:
        atom = pair[0]
        if atom not in anion:
            at1.append(atom)
        else:
            at2.append(atom)
    return list(set(at1)), list(set(at2))


def main():
    func = 'St'
    if func == 'Ct':
        top = "../md.tpr"
        uta = mda.Universe(top)
        aP = uta.select_atoms('type opls_480')  # fetch anions
        lP = np.linspace(0, len(aP) - 1, len(aP), dtype = int)
        ac = h5py.File('../h5ac/pils.h5', 'r')
        al = h5py.File('../h5al/pils.h5', 'r')
        h5 = h5py.File('aniongroups.h5', 'w')
        dtype = h5py.special_dtype(vlen = np.dtype('int32'))
        ds_AT1 = get_or_create_dataset(h5, 'AT1', dtype)
        ds_AT2 = get_or_create_dataset(h5, 'AT2', dtype)
        ds_AT3 = get_or_create_dataset(h5, 'AT3', dtype)
        ASSOac = rawLoad(ac, 'htt')
        ASSOal = rawLoad(al, 'htt')
        frames = 10001
        for time in range(frames):
            if time % 100 == 0:
                print("Processing {:10.3f}".format((time + 1) * 1.0 / frames))
            polyIL, Both = getgroups(time, ASSOac, ASSOal)
            Vehicle = np.setdiff1d(lP, polyIL + Both)
            ds_append(ds_AT1, polyIL)
            ds_append(ds_AT2, Both)
            ds_append(ds_AT3, Vehicle)
        h5.close()
        h5file = 'aniongroups.h5'
        h5tag = 'AT1'
        anionACF = finalGetC_i(h5file, h5tag, 0, 1, 48)
        anionPrint('at1ACF', anionACF)
        h5tag = 'AT2'
        anionACF = finalGetC_i(h5file, h5tag, 0, 1, 48)
        anionPrint('at2ACF', anionACF)
        h5tag = 'AT3'
        anionACF = finalGetC_i(h5file, h5tag, 0, 1, 48)
        anionPrint('at3ACF', anionACF)
    elif func == 'St':
        top = "../md.tpr"
        uta = mda.Universe(top)
        aP = uta.select_atoms('type opls_480')  # fetch anions
        lP = np.linspace(0, len(aP) - 1, len(aP), dtype = int)
        ac = h5py.File('../h5ac/pils.h5', 'r')
        al = h5py.File('../h5al/pils.h5', 'r')
        h5 = h5py.File('aniongroups.h5', 'w')
        dtype = h5py.special_dtype(vlen = np.dtype('int32'))
        ds_AT1 = get_or_create_dataset(h5, 'AT1', dtype)
        ds_AT2 = get_or_create_dataset(h5, 'AT2', dtype)
        ds_AT3 = get_or_create_dataset(h5, 'AT3', dtype)
        ASSOac = rawLoad(ac, 'htt')
        ASSOal = rawLoad(al, 'htt')
        frames = 10001
        for time in range(frames):
            if time % 100 == 0:
                print("Processing {:10.3f}".format((time + 1) * 1.0 / frames))
            polyIL, Both = getgroups(time, ASSOac, ASSOal)
            Vehicle = np.setdiff1d(lP, polyIL + Both)
            ds_append(ds_AT1, polyIL)
            ds_append(ds_AT2, Both)
            ds_append(ds_AT3, Vehicle)
            if time == 2000:
                break
        h5.close()
        h5file = 'aniongroups.h5'
        h5tag = 'AT1'
        anionACF = finalGetC_c(h5file, h5tag, 0, 1, 48)
        anionPrint('at1ACFSt', anionACF)
        h5tag = 'AT2'
        anionACF = finalGetC_c(h5file, h5tag, 0, 1, 48)
        anionPrint('at2ACFSt', anionACF)
        h5tag = 'AT3'
        anionACF = finalGetC_c(h5file, h5tag, 0, 1, 48)
        anionPrint('at3ACFSt', anionACF)        


if __name__ == "__main__":
    main()
