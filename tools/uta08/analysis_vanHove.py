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


def alpha2Print(filename, alpha2):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# alpha2', file=anaout)
        for i in range(len(alpha2)):
            print('{:10.5f}'.format(alpha2[i]), file=anaout)


def ds_append(ds, val):
    ds.resize(ds.shape[0]+1, axis=0)
    ds[-1] = val


def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def getGself(COORDs, t0, t):
    """
    The dynamical heterogeneity deviates from the ideal Gaussian parameter
    J. Phys. Chem. B 2021, 125, 372−381
    """
    gs = []
    r = np.linspace(0, 20, 1001)
    displ = np.sqrt(np.sum(np.square(COORDs[t] - COORDs[t0]), axis = 1))
    Natoms = len(displ)
    for ri in r:
        temp = [1.0 for rj in displ if np.abs(ri - rj) <= 1e-3]
        gs.append(np.sum(temp) / Natoms)   
    return gs


def Gself(COORDs, tstar, tf, reservoir):
    a = 0
    count = 0
    for t0 in reservoir:
        if t0 + tstar <= tf:
            a += getGself(COORDs, t0, t0 + tstar)
            count += 1
        else:
            break
    return a / count


def finalGself(COORDs, tstar, trestart, nt):
    """
    This function gets the final data of the datapoint_i graph.
    """
    tf = len(COORDs) - 1
    num_points = len(COORDs)
    pool = mp.Pool(processes=nt)
    reservoir = [x * trestart for x in range(tf) if x * trestart < tf]
    func = functools.partial(Gself, COORDs, tstar, tf)
    return_data = pool.map(func, reservoir)
    return return_data


def main():
    # generic parameters
    nframes = 10001                        # number of frames
    nt = 56                                # number of processors
    trestart = 10                          # frequency for picking up the reference point
    tstar = 1500                           # t star obtained from the alpha2 parameter
    # generate the coordinates object
    top = "../cg_topol.tpr"
    trj = "unwrap.xtc"
    atp = "PF"                             # fetch central particle
    uta = mda.Universe(top, trj)
    atg = uta.select_atoms("type " + atp)  # atom group
    COORDs = []
    for ts in uta.trajectory:
        if ts.frame >= 0:                  # start from arbitrary given reference frame
            COORDs.append(atg.positions)
        if ts.frame % 100 == 0:
            print("currently processing the frame {}".format(ts.frame))
        if ts.frame == nframes:
            break
    alpha2 = finalGself(COORDs, tstar, trestart, nt)
    alpha2Print('alpha2', alpha2)

if __name__ == "__main__":
    main()
