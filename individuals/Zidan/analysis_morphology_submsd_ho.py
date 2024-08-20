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


def msdPrint(filename, time, msd0):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# time total', file=anaout)
        for i in range(len(time)):
            print('{:10.5f} {:10.5f}'.format(time[i], msd0[i]), file=anaout)


def timestamp(TIME):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(TIME + " =", current_time)


def main():
    timestamp("Start Time")
    # load trajectory for MSD
    top = "cg_topol.tpr" 
    trj = "unwrap.xtc" 
    uta = mda.Universe(top, trj) 
    aP = uta.select_atoms('type PF')  # fetch anions
    START = np.linspace(0, 8000, 101, dtype = int)
    tau = 2000
    MSDg0 = []
    for start in START:
        print(start)
        COORDs = []
        for ts in uta.trajectory:
            if ts.frame >= start:
                COORDs.append(aP.positions)
            if ts.frame == start + tau + 1:
                break
        dispg0 = np.zeros((tau + 1, len(aP)))
        frame = 0
        for time in range(start, tau + start + 1):
            displ = COORDs[frame] - COORDs[0]
            msd = np.sum(np.square(displ), axis = 1)
            dispg0[frame] = msd
            frame += 1
        msdg0 = np.zeros(tau + 1)
        for time in range(1, tau + 1):
            msdg0[time] = np.average(dispg0[time])
        MSDg0.append(msdg0)
    msd0 = np.zeros(tau + 1)
    for time in range(1, tau + 1):
        tmp0 = []
        for i in range(len(START)):
            tmp0.append(MSDg0[i][time])
        tmp0 = np.asarray(tmp0)
        msd0[time] = np.average(tmp0)
    xdata = np.linspace(0, tau, tau + 1)
    msdPrint('subMSD', xdata, msd0)
    timestamp("End Time")
    print("Sub-MSD analysis finished!")

if __name__ == "__main__":
    main()
