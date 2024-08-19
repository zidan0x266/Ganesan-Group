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


def msdPrint(filename, time, msdg1, msdg2):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# time msdg1 msdg2', file=anaout)
        for i in range(len(time)):
            print('{:10.5f} {:10.5f} {:10.5f}'.format(time[i], msdg1[i], msdg2[i]), file=anaout)


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
    START = np.linspace(0, 9600, 1001, dtype = int)
    tau = 200
    # load ion pair associations
    ac = h5py.File('../test/ac.h5', 'r')
    al = h5py.File('../test/al.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')   
    MSDg1, MSDg2 = [], []
    for start in START:
        print(start)
        COORDs = []
        for ts in uta.trajectory:
            if ts.frame >= start:
                COORDs.append(aP.positions)
            if ts.frame == start + tau + 1:
                break
        dispg1, dispg2 = np.zeros((tau + 1, len(aP))), np.zeros((tau + 1, len(aP)))
        stateg1 = np.zeros(len(aP))
        stateg2 = np.zeros(len(aP))
        firstframe = True
        frame = 0
        for time in range(start, tau + start + 1):
            atomg1, atomg2 = getgroups(time, ASSOac, ASSOal)[0:2]
            #atomg3 = np.setdiff1d(lP, atomg1 + atomg2)  # ignore this contribution
            if firstframe:  # marker the inital status of the ion group
                for atom in atomg1:
                    stateg1[atom] = 1
                for atom in atomg2:
                    stateg2[atom] = 1
                firstframe = False
            tmpg1 = np.zeros(len(aP))
            tmpg2 = np.zeros(len(aP))
            for atom in atomg1:   
                tmpg1[atom] = 1
            for atom in atomg2:
                tmpg2[atom] = 1
            stateg1 *= tmpg1
            stateg2 *= tmpg2
            displ = COORDs[frame] - COORDs[0]
            msd = np.sum(np.square(displ), axis = 1)
            dispg1[frame] = msd * stateg1
            dispg2[frame] = msd * stateg2
            frame += 1
        msdg1 = np.zeros(tau + 1)
        msdg2 = np.zeros(tau + 1)
        for time in range(1, tau + 1):
            msdg1[time] = dispg1[time][np.nonzero(dispg1[time])].mean()
            msdg2[time] = dispg2[time][np.nonzero(dispg2[time])].mean()
        MSDg1.append(msdg1)
        MSDg2.append(msdg2)
    msd1 = np.average(MSDg1, axis = 0)
    msd2 = np.average(MSDg2, axis = 0)
    xdata = np.linspace(0, tau, tau + 1)
    msdPrint('subGroupMSD', xdata, msd1, msd2)

if __name__ == "__main__":
    main()