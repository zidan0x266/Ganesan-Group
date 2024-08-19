"""
Copyright (C) 2018-2022 Zidan Zhang <zhangzidan@gmail.com>

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


def displPrint(filename, displacement, atom):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# X Y Z', file=anaout)
        for i in range(len(displacement)):
            print('{:10.5f} {:10.5f} {:10.5f}'.format(displacement[i][atom][0], displacement[i][atom][1], displacement[i][atom][2]), file=anaout)


def ds_append(ds, val):
    ds.resize(ds.shape[0]+1, axis=0)
    ds[-1] = val


def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def getdp_i(COORDs, t0, t):
    """
    The dynamical heterogeneity deviates from the ideal Gaussian parameter
    J. Phys. Chem. B 2021, 125, 372−381
    """
    displ = COORDs[t] - COORDs[t0]
    return displ


def main():
    # generic parameters
    nframes = 200001                         # number of frames
    natoms = 256
    # generate the coordinates object
    top = "cg_topol.tpr"
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
    displ = []
    for time in range(nframes):
        displ.append(getdp_i(COORDs, 0, time))
    for atom_id in range(natoms):
        if atom_id % 100 == 0:
            print("currently processing the atom {}".format(atom_id))
        displPrint('pf6_' + str(atom_id), displ, atom_id)


if __name__ == "__main__":
    main()
