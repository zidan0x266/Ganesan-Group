"""
Copyright (C) 2018-2021 Zidan Zhang <zhangzidan@gmail.com>

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
import numpy as np
import h5py
import time
import MDAnalysis as mda

def construction(h5, dataset_name, length, dtype):
    dtype = h5py.vlen_dtype(np.dtype('int32'))
    if dataset_name not in h5:
        #return h5.create_dataset(dataset_name, (0, length), dtype=dtype)
        return h5.create_dataset(dataset_name, (0, length), maxshape=(None, length), dtype=dtype, chunks=(100, length))
    return h5[dataset_name]


def ds_append(ds, val):
    ds.resize(ds.shape[0]+1, axis=0)
    ds[-1] = val


def rawSave(h5, dataset_name, raw_data, length):
    dtype = h5py.vlen_dtype(np.dtype('int32'))
    if dataset_name not in h5:
        dataset = h5.create_dataset(dataset_name, (len(raw_data), length), dtype=dtype)
    else:
        dataset = h5[dataset_name]
    for i, v in enumerate(raw_data):
        dataset[i] = v


def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def acfCtPrint(filename, aacf):
    """
    This function print the intermittent Association Lifetime
    """
    acfcout = open(filename + '_Ct.dat', 'w')
    print('# Time X_t', file=acfcout)
    tc = np.linspace(0, len(aacf) - 1, len(aacf))
    for i in range(0, len(aacf)):
        print('{} {:10.5f}'.format(tc[i], aacf[i]), file=acfcout)
    acfcout.close()


def getC_i(Vehicle, t0, t):
    """
    This function give the intermitent Association Lifetime
    C_i = <h(t0)h(t)>/<h(t0)> between t0 and t
    """
    if Vehicle[t0].shape[0] == 0:
        return 0.0

    origin = [len(x) for x in Vehicle[t0]]
    overlap = [len(np.intersect1d(x, y)) for x, y in zip(Vehicle[t0], Vehicle[t])]
    temp = [x * 1.0 / y for x, y in zip(overlap, origin) if y != 0]
    C_i = np.average(temp)

    return C_i


def intervC_i(h5filename, dataset_name, tf, reservoir, dt):
    """
    This function gets all the data for the h(t0)h(t0+dt), where
    t0 = 1,2,3,...,tf. This function give us a point of the final plot
    C(t) vs t,
    """
    Vehicle, h5 = rawLoad(h5filename, dataset_name)
    a = 0
    count = 0
    for t0 in reservoir:
        if t0 + dt <= tf:
            a += getC_i(Vehicle, t0, t0 + dt)
            count += 1
    h5.close()
    return a / count


def finalGetC_i(h5filename, dataset_name, trestart, nt):
    """
    This function gets the final data of the C_i graph.
    """
    Vehicle, h5 = rawLoad(h5filename, dataset_name)
    tf = Vehicle.shape[0] - 1
    num_Vehicle = Vehicle.shape[0]
    h5.close()

    pool = mp.Pool(processes=nt)
    reservoir = [x * trestart for x in range(tf) if x * trestart < tf]
    func = functools.partial(intervC_i, h5filename, dataset_name, tf, reservoir)
    return_data = pool.map(func, list(range(num_Vehicle)))
    return return_data


def get_solvation(cur_atom, ASSO):
    """
    get the first solvation shell for each lithium of the single frame
    """
    atom = cur_atom
    sol_shell = []  # first solvation shell
    for pair in range(len(ASSO)):
        atom1, atom2 = ASSO[pair]
        if atom1 == atom:
            sol_shell.append(atom2)
    return list(set(sol_shell))


def solvation(pst, natoms, frame_id):
    """
    get the first solvation shell for each lithium
    """
    H5file = h5py.File(pst, 'r')
    ASSO = rawLoad(H5file, 'htt')
    Sol_shell = []
    for cur_atom in range(natoms):
        sol_shell = get_solvation(cur_atom, ASSO[frame_id])
        Sol_shell.append(sol_shell)
    return Sol_shell


def main():
    func = 0
    frames = 2001  # number of frames
    nt = 56  # number of processors
    trestart = 10  # frequency of the reference point
    top = "../cg_topol.tpr"
    pst = "../asso/pils.h5"
    atp = "PF"  # fetch central particle
    uta = mda.Universe(top)
    atg = uta.select_atoms("type " + atp)
    if func == 0:
        dtype = h5py.vlen_dtype(np.dtype('int32'))
        h5 = h5py.File('solvation.h5', 'w')  # save vehicle information
        dataset = construction(h5, 'shell', len(atg), dtype)
        frame_ids = [frame for frame in range(frames)]
        solvation_analyse = functools.partial(solvation, pst, len(atg))
        pool = mp.Pool(processes=nt)
        data = pool.imap(solvation_analyse, frame_ids)
        for Solvation in data:
            ds_append(dataset, Solvation)
        h5.close()
    C_t = finalGetC_i('solvation.h5', 'shell', trestart, nt)
    acfCtPrint('Solvation', C_t)
    

if __name__ == "__main__":
    main()
