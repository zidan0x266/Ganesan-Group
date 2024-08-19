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

import h5py
import numpy as np
import MDAnalysis as mda

def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def msdPrint(filename, time, msd):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# time msd', file=anaout)
        for i in range(len(time)):
            print('{:10.5f} {:10.5f}'.format(time[i], msd[i]), file=anaout)


def getgroups(frame, ASSOac, ASSOal):  # get the co-coordination anions
    type1, type2 = [], []
    at1, at2 = [], []
    anion = list(set([x[0] for x in ASSOal[frame]]))  # find the anions in co-coordination
    for pair in ASSOac[frame]:
        atom = pair[0]
        if atom not in anion:
            type1.append(pair)  # anions only associated with cation
            at1.append(atom)
        else:
            type2.append(pair)  # anions in the co-coordination
            at2.append(atom)
    return list(set(at1)), list(set(at2)), type1, type2


def main():
    # load trajectory for MSD
    top = "md.tpr" 
    trj = "unwrap.xtc" 
    uta = mda.Universe(top, trj) 
    aP = uta.select_atoms('type opls_480')  # fetch anions
    #lP = np.linspace(0, len(aP) - 1, len(aP), dtype = int)  # create the atom index
    START = np.linspace(0, 9000, 10, dtype = int)
    tau = 1500
    # load ion pair associations
    ac = h5py.File('../test/ac.h5', 'r')
    al = h5py.File('../test/al.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')   
    MSDs = []
    for start in START:
        print(start)
        COORDs = []
        for ts in uta.trajectory:
            if ts.frame >= start:
                COORDs.append(aP.positions)
            if ts.frame == start + tau + 1:
                break
        displ = np.zeros((tau + 1, len(aP)))
        state = np.zeros(len(aP))
        msds = np.zeros(tau + 1)
        firstframe = True
        frame = 0
        for time in range(start, tau + start + 1):
            anions = getgroups(time, ASSOac, ASSOal)[0]
            #atomg3 = np.setdiff1d(lP, atomg1 + atomg2)  # ignore this contribution
            if firstframe:  # marker the inital status of the ion group
                for atom in anions:
                    state[atom] = 1
                firstframe = False
            tmp = np.zeros(len(aP))
            for atom in anions:   
                tmp[atom] = 1
            state *= tmp
            disr = COORDs[frame] - COORDs[0]
            msdt = np.sum(np.square(disr), axis = 1)
            displ[frame] = msdt * state
            if frame != 0:
                msds[frame] = displ[frame][np.nonzero(displ[frame])].mean()
            frame += 1
        MSDs.append(msds)
    msdf = np.average(MSDs, axis = 0)
    xdata = np.linspace(0, tau, tau + 1)
    msdPrint('anionMSD', xdata, msdf)

if __name__ == "__main__":
    main()