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


def travelPrint(filename, distance):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# distance', file=anaout)
        for i in range(len(distance)):
            print('{:10.5f}'.format(distance[i]), file=anaout)


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
        #return h5.create_dataset(dataset_name, (0, length), dtype=dtype)
        return h5.create_dataset(dataset_name, (0, length), maxshape=(None, length), dtype=dtype, chunks=(100, length))
    return h5[dataset_name]


def get_solvation(lithium, ASSOal):
    """
    get the first solvation shell for each lithium of the single frame
    """
    atom = lithium
    sol_shell = []  # first solvation shell
    for pair in range(len(ASSOal)):
        atom1, atom2 = ASSOal[pair]
        if atom2 == atom:
            sol_shell.append(atom1)
    return list(set(sol_shell))


def solvation(pst, Nlithiums, frame_id):
    """
    get the first solvation shell for each lithium
    """
    H5al = h5py.File(pst, 'r')
    ASSOal = rawLoad(H5al, 'htt')
    Sol_shell = []
    for mylithium in range(Nlithiums):
        sol_shell = get_solvation(mylithium, ASSOal[frame_id])
        Sol_shell.append(sol_shell)
    return Sol_shell


def generate_solvation(filename, pst, Nlithiums, nframes, nt):
    dtype = h5py.vlen_dtype(np.dtype('int32'))
    h5 = h5py.File(filename + '.h5', 'w')  # save solvation shell information
    Sol_shell = construction(h5, 'shell', Nlithiums, dtype)
    frame_ids = [frame for frame in range(nframes)]
    solvation_analyse = functools.partial(solvation, pst, Nlithiums)
    pool = mp.Pool(processes=nt)
    data = pool.imap(solvation_analyse, frame_ids)
    for Shell in data:
        ds_append(Sol_shell, Shell)
    h5.close()


def get_prob(old, new):
    if list(new) == list(old):
        prob = 1.0
    else:
        prob = 0.0
    return prob


def main():
    # generic parameters
    func = 0  # if 0: generate the .h5 file, else using exsiting .h5 file
    nframes = 10001  # number of frames
    nt = 56  # number of processors
    # generate the coordinates object
    top = "../md.tpr"
    trj = "../subMSDs/unwrap.xtc"
    pst = "../h5al/pils.h5"
    atp = "opls_406"  # fetch lithium atoms
    uta = mda.Universe(top, trj)
    atg = uta.select_atoms("type " + atp)
    #lL = np.linspace(0, len(atg) - 1, len(atg), dtype = int)  # create the atom index
    if func == 0:
        generate_solvation('solvation', pst, len(atg), nframes, nt)
    START = np.linspace(0, 5000, 6, dtype = int)
    tau = 500  # number of frames for analysis the cumulative travel distance
    Vehicle, h5 = rawLoad('solvation.h5', 'shell')
    # load ion pair associations
    cumulative_length = []
    for start in START:
        print("The current reference frame is: " + str(start))
        vehicle_state = np.zeros((len(atg), tau + 1))  # initialize the lithium solvation shell state function
        initialization = True
        counter = 0
        for time in range(start, tau + start + 1):
            if initialization:
                old_sol_shell = Vehicle[start]
                new_sol_shell = Vehicle[start]
                v_state = [get_prob(old, new) for old, new in zip(old_sol_shell, new_sol_shell)]
                vehicle_state[:][:, counter] = v_state
                counter += 1
                initialization = False
            else:
                new_sol_shell = Vehicle[time]
                v_state = [get_prob(old, new) for old, new in zip(old_sol_shell, new_sol_shell)]
                for i in range(len(v_state)):
                    if vehicle_state[i][counter - 1] == 1:
                        vehicle_state[i][counter] = v_state[i]
                    else:
                        continue
                old_sol_shell = new_sol_shell
                counter += 1
        COORDs = []
        for ts in uta.trajectory:
            if ts.frame >= start:
                COORDs.append(atg.positions)
            if ts.frame == start + tau + 1:
                break
        displacement = np.zeros((len(atg), tau + 1))
        accumulation = np.zeros(len(atg))
        frame = 0
        for time in range(start, tau + start + 1):
            if frame == 0 or frame == 1:
                displ = COORDs[frame] - COORDs[0]
                squared = np.sum(np.square(displ), axis = 1)
                accumulation += np.sqrt(squared)
                displacement[:][:, frame] = accumulation
            else:
                displ = np.abs(COORDs[frame] - COORDs[frame - 1])
                squared = np.sum(np.square(displ), axis = 1)
                accumulation += np.sqrt(squared)
                displacement[:][:, frame] = accumulation * list(vehicle_state[:][:, frame])
            frame += 1
        for atom in range(len(atg)):
            distance = [dis for dis in displacement[atom] if dis > 0]
            if len(distance) != 0:
                cumulative_length.append(distance[-1])
            else:
                cumulative_length.append(0.0)
    travelPrint('cumulative_length', cumulative_length)

if __name__ == "__main__":
    main()
