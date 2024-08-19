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
import h5py

import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
import MDAnalysis.lib.distances as distances
import numpy as np
from collections import Counter

"""
This program is to do the hopping analysis that was proposed in
JPCP 2010, 114, 877-881.
"""

def ds_append(ds, val):
    ds.resize(ds.shape[0]+1, axis=0)
    ds[-1] = val


def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def construction(h5, dataset_name, length, dtype):
    dtype = h5py.vlen_dtype(np.dtype('int32'))
    if dataset_name not in h5:
        return h5.create_dataset(dataset_name, (0, length), maxshape=(None, length), dtype=dtype, chunks=(100, length))
    return h5[dataset_name]


def solvation(top, trj, ag1, ag2, cutoff, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ag1)
    group2 = uta.select_atoms("type " + ag2)
    ts = uta.trajectory[frame_id]
    cell = ts.dimensions
    asso = []
    for atom in range(len(group1)):
        distpair = distances.capped_distance(group1.positions[atom], group2.positions, cutoff, box = cell)[0]
        shell = [pair[1] for pair in distpair]
        asso.append(shell)
    return asso


def generate_solvation(filename, natoms, frame_ids, top, trj, ag1, ag2, cutoff, nt):
    dtype = h5py.vlen_dtype(np.dtype('int32'))
    h5 = h5py.File(filename + '.h5', 'w')  # save solvation shell information
    Sol_shell = construction(h5, 'shell', natoms, dtype)
    #frame_ids = [frame for frame in range(nframes)]
    solvation_analyse = functools.partial(solvation, top, trj, ag1, ag2, cutoff)
    pool = mp.Pool(processes=nt)
    data = pool.imap(solvation_analyse, frame_ids)
    for Shell in data:
        ds_append(Sol_shell, Shell)
    h5.close()


def get_delta(pre, old, new):
    if list(new) == list(old):
        delta = 0.0
    elif list(new) != list(old) and list(new) == list(pre):
        delta = -1.0
    else:
        delta = 1.0
    return delta


def main():
    load_traj = True
    create_solvation = False
    h5filename = "state_"
    nt=20  # number of processors
    if load_traj:
        top = "../gentpr/cg_topol.tpr"
        trj = "../data_100ns.xtc"
        atom_group1 = "TF"
        atom_group2 = "IM"
        atom_group3 = "LI"
        cutoff1 = 7.7  # cutoff for anion - cation
        cutoff2 = 5.1  # cutoff for anion - lithium
        uta = mda.Universe(top, trj)
        natoms = len(uta.select_atoms("type " + atom_group3))  # number of atoms in the center of the solvation shell
        frame_ids = [ts.frame for ts in uta.trajectory]
    if create_solvation:  # by default create the .h5 file that contains the solvation shell
        generate_solvation('solvation', natoms, frame_ids, top, trj, atom_group3, atom_group1, cutoff2, nt)
    Vehicle, h5 = rawLoad('solvation.h5', 'shell')  # load saved vehicle information
    tau = 5  # number of frames for analysis the cumulative travel distance
    length = int((len(Vehicle) - 1) / tau) + 1
    ht = np.zeros((length, natoms))
    firstframe, secondframe = True, False
    for i in range(length):
        if firstframe:
            pre_sol_shell = Vehicle[i * tau]
            old_sol_shell = Vehicle[i * tau]
            new_sol_shell = Vehicle[i * tau]
            v_state = [get_delta(pre, old, new) for pre, old, new in zip(pre_sol_shell, old_sol_shell, new_sol_shell)]
            ht[i] = v_state
            firstframe = False
            secondframe = True
        elif secondframe:
            pre_sol_shell = Vehicle[(i - 1)* tau]
            old_sol_shell = Vehicle[(i - 1)* tau]
            new_sol_shell = Vehicle[i * tau]
            v_state = [get_delta(pre, old, new) for pre, old, new in zip(pre_sol_shell, old_sol_shell, new_sol_shell)]
            ht[i] = ht[i - 1] + v_state
            secondframe = False
        else:
            pre_sol_shell = Vehicle[(i - 2)* tau]
            old_sol_shell = Vehicle[(i - 1)* tau]
            new_sol_shell = Vehicle[i * tau]
            v_state = [get_delta(pre, old, new) for pre, old, new in zip(pre_sol_shell, old_sol_shell, new_sol_shell)]
            ht[i] = ht[i - 1] + v_state
    statef = h5py.File(h5filename + str(tau) + ".h5", "w")
    dataset = statef.create_dataset("aht", (length, natoms), dtype = 'f', data = ht)  # atom specific h(t)
    statef.close()


if __name__ == "__main__":
    main()