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
from datetime import datetime

def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def msdPrint(filename, time, msd0, msd1, msd2):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# time total bulk interface', file=anaout)
        for i in range(len(time)):
            print('{:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(time[i], msd0[i], msd1[i], msd2[i]), file=anaout)


def timestamp(TIME):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(TIME + " =", current_time)


def getgroups(frame, ASSObk, ASSOif):  # get the co-coordination anions
    abk = list(set([x[0] for x in ASSObk[frame]]))
    aif = list(set([x[0] for x in ASSOif[frame]]))
    return list(set(abk)), list(set(aif))


def main():
    timestamp("Start Time")
    # load trajectory for MSD
    top = "cg_topol.tpr" 
    trj = "unwrap.xtc" 
    uta = mda.Universe(top, trj) 
    aP = uta.select_atoms('type PF')  # fetch anions
    lP = np.linspace(0, len(aP) - 1, len(aP), dtype = int)  # create the atom index
    START = np.linspace(0, 8000, 101, dtype = int)
    tau = 2000
    # load ion pair associations
    h5 = h5py.File('../../dynamics/la01/pils.h5', 'r')
    ASSObk = rawLoad(h5, 'in_htt')
    ASSOif = rawLoad(h5, 'it_htt')   
    MSDg0, MSDg1, MSDg2 = [], [], []
    for start in START:
        print(start)
        COORDs = []
        for ts in uta.trajectory:
            if ts.frame >= start:
                COORDs.append(aP.positions)
            if ts.frame == start + tau + 1:
                break
        dispg0 = np.zeros((tau + 1, len(aP)))
        dispg1 = np.zeros((tau + 1, len(aP)))
        dispg2 = np.zeros((tau + 1, len(aP)))
        stateg1 = np.zeros(len(aP))
        stateg2 = np.zeros(len(aP))
        firstframe = True
        frame = 0
        for time in range(start, tau + start + 1):
            if firstframe:  # marker the inital status of the ion group
                atomg1, atomg2 = getgroups(time, ASSObk, ASSOif)
                for atom in atomg1:
                    stateg1[atom] = 1
                for atom in atomg2:
                    stateg2[atom] = 1
                firstframe = False
            displ = COORDs[frame] - COORDs[0]
            msd = np.sum(np.square(displ), axis = 1)
            dispg0[frame] = msd
            dispg1[frame] = msd * stateg1
            dispg2[frame] = msd * stateg2
            frame += 1
        msdg0 = np.zeros(tau + 1)
        msdg1 = np.zeros(tau + 1)
        msdg2 = np.zeros(tau + 1)
        for time in range(1, tau + 1):
            msdg0[time] = np.average(dispg0[time])
            msdg1[time] = dispg1[time][np.nonzero(dispg1[time])].mean()
            msdg2[time] = dispg2[time][np.nonzero(dispg2[time])].mean()
        MSDg0.append(msdg0)
        MSDg1.append(np.nan_to_num(msdg1))
        MSDg2.append(np.nan_to_num(msdg2))
    msd0 = np.zeros(tau + 1)
    msd1 = np.zeros(tau + 1)
    msd2 = np.zeros(tau + 1)
    for time in range(1, tau + 1):
        tmp0, tmp1, tmp2 = [], [], []
        for i in range(len(START)):
            tmp0.append(MSDg0[i][time])
            tmp1.append(MSDg1[i][time])
            tmp2.append(MSDg2[i][time])
        tmp0 = np.asarray(tmp0)
        tmp1 = np.asarray(tmp1)
        tmp2 = np.asarray(tmp2)
        msd0[time] = np.average(tmp0)
        msd1[time] = tmp1[np.nonzero(tmp1)].mean()
        msd2[time] = tmp2[np.nonzero(tmp2)].mean()
    xdata = np.linspace(0, tau, tau + 1)
    msdPrint('subMSD', xdata, msd0, msd1, msd2)
    timestamp("End Time")
    print("Sub-MSD analysis finished!")

if __name__ == "__main__":
    main()
