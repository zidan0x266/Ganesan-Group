"""
Copyright (C) 2018-2023 Zidan Zhang <zhangzidan@gmail.com>

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
from collections import Counter


def assoFreq(frequence):
    """
    This function gets the normolized probability for association events.
    """
    fNum = frequence[0]
    for i in range(1, len(frequence)):
        fNum += frequence[i]
    numAsso = {}
    for k, v in fNum.items():
        numAsso[k] = v / float(len(frequence))
    totalAsso = sum(numAsso.values())
    avgAsso = {}
    for k, v in numAsso.items():
        avgAsso[k] = v / float(totalAsso)
    pairNum = []
    pairPro = []
    for k in sorted(avgAsso.keys()):
        pairNum.append(k)
        pairPro.append(avgAsso[k])
    return pairNum, pairPro


def anaPrint(filename, arrayNum, arrayPro):  # print association properties
    anaout = open(filename + '.dat', 'w')
    print('# ' + filename + ' Probability', file=anaout)
    for i in range(0, len(arrayNum)):
        print('{} {:10.5f}'.format(arrayNum[i], arrayPro[i]), file=anaout)
    anaout.close()


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


def get_solvation(lithium, ASSO):
    """
    get the first solvation shell for each lithium of the single frame
    """
    atom = lithium
    sol_shell = []  # first solvation shell
    for pair in range(len(ASSO)):
        atom1, atom2 = ASSO[pair]
        if atom1 == atom:
            sol_shell.append(atom2)
    return list(set(sol_shell))


def solvation(pst, Nlithiums, frame_id):
    """
    get the first solvation shell for each lithium
    """
    H5asso = h5py.File(pst, 'r')
    ASSO = rawLoad(H5asso, 'htt')
    Sol_shell = []
    for mylithium in range(Nlithiums):
        sol_shell = get_solvation(mylithium, ASSO[frame_id])
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


def get_prob(old, new, probability):
    differ = np.intersect1d(new, old)
    ratio = float(len(differ)) / len(old)
    if ratio >= probability:
        prob = 1.0
    else:
        prob = 0.0
    return prob


def getC_c(ASP, t0, t):
    """
    This function give the continous Association Lifetime
    C_c = <h(t0)H(t)>/<h(t0)> between t0 and t.
    """
    if t0 == t:
        return 1.0

    dt = 1
    ht0 = ASP[t0]
    htt_old = ASP[t0]

    while t0 + dt <= t:
        htt = ASP[t0 + dt]
        Htt = np.intersect1d(htt_old, htt)
        # Htt = set(htt_old) & set(htt)  # comparison between current frame and previous frame
        t0 += dt
        htt_old = Htt
    C_c = len(Htt) / float(len(ht0))
    if ASP[t0].shape[0] == 0:
        return 0
    else:
        return C_c


def intervC_c(h5filename, dataset_name, t0, tf, reservoir, dt):
    """
    This function gets all the data for the h(t0)H(t0+dt), where
    t0 = 1,2,3,...,tf. This function give us one point of the final plot
    S(t) vs t.
    """
    ASP, h5 = rawLoad(h5filename, dataset_name)
    a = 0
    count = 0
    for t0 in reservoir:
        if t0 + dt <= tf:
            if t0 == t0 + dt:
                break
            a += getC_c(ASP, t0, t0 + dt)
            count += 1

    h5.close()

    if count == 0:
        return 1.0

    return a / count


def finalGetC_c(h5filename, dataset_name, t0, trestart, nt):
    """
    This function gets the final data of the C_c graph.
    """
    ASP, h5 = rawLoad(h5filename, dataset_name)
    tf = ASP.shape[0] - 1
    num_ASP = ASP.shape[0]
    print(('finalGetC_c tf: {}'.format(tf)))
    h5.close()

    pool = mp.Pool(processes=nt)
    reservoir = [x * trestart for x in range(tf) if x * trestart < tf]
    func = functools.partial(intervC_c, h5filename, dataset_name, t0, tf, reservoir)
    return_data = pool.map(func, list(range(num_ASP)))
    return return_data


def main():
    # generic parameters
    func = 1  # if 0: generate the .h5 file, else using exsiting .h5 file
    nframes = 10001  # number of frames
    nt = 56  # number of processors
    # generate the coordinates object
    top = "../cg_topol.tpr"
    trj = "unwrapped.xtc"
    pst = "../h5la/pils.h5"
    atp = "LI"  # fetch central particle
    uta = mda.Universe(top, trj)
    atg = uta.select_atoms("type " + atp)
    #lL = np.linspace(0, len(atg) - 1, len(atg), dtype = int)  # create the atom index
    if func == 0:
        generate_solvation('solvation', pst, len(atg), nframes, nt)
    tau = 4500  # number of frames for analysis the cumulative travel distance
    Vehicle, h5 = rawLoad('../h5la/solvation.h5', 'shell')
    # Get coordinates
    COORDs = []
    freqVehicle = []
    for ts in uta.trajectory:
        COORDs.append(atg.positions)
    START = np.linspace(0, 9000, 91, dtype = int)
    for start in START:
        # Here the displacements of individual atom will be calculated without statistical enhancement
        atoms = np.arange(len(atg))
        atoms2rm = []
        remainassos = []
        frames = np.linspace(start, 10000, 10001 - start, dtype = int)
        for frame in frames:
            displ = COORDs[frame] - COORDs[start]
            msd = np.sum(np.square(displ), axis = 1)
            displacement = np.sqrt(msd)      
            atoms = np.setdiff1d(atoms, atoms2rm)
            if len(atoms) == 0:
                break
            else:
                for atom in atoms:
                    distance = displacement[atom]
                    if distance >= 5.0:
                        atoms2rm.append(atom)
                        remains = np.intersect1d(Vehicle[frame][atom], Vehicle[start][atom])
                        remainassos.append(len(remains))
        print(remainassos)
        freqVehicle.append(Counter(remainassos))
    Num, Pro = assoFreq(freqVehicle)
    anaPrint('remainPair', Num, Pro)

if __name__ == "__main__":
    main()