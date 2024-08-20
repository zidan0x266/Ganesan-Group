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
    START = np.linspace(0, 9700, 1000, dtype = int)
    tau = 200  
    MSDs = []
    for start in START:
        print(start)
        COORDs = []
        for ts in uta.trajectory:
            if ts.frame >= start:
                COORDs.append(aP.positions)
            if ts.frame == start + tau + 1:
                break
        msds = np.zeros(tau + 1)
        firstframe = True
        frame = 0
        for time in range(start, tau + start + 1):
            displ = COORDs[frame] - COORDs[0]
            tmp = np.sum(np.square(displ), axis = 1)
            if frame != 0:
                msds[frame] = np.average(tmp)
            frame += 1
        MSDs.append(msds)
    msdfinal = np.average(MSDs, axis = 0)
    xdata = np.linspace(0, tau, tau + 1)
    msdPrint('GroupMSD', xdata, msdfinal)

if __name__ == "__main__":
    main()