"""
Copyright (C) 2018-2024 Zidan Zhang <zhangzidan@gmail.com>

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


def msdPrint(filename, time, msd0, msd1_to, msd1_bk, msd1_if):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# time, total, msd1_to, msd1_bk, msd1_if', file=anaout)
        for i in range(len(time)):
            print('{:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(time[i], msd0[i], msd1_to[i], msd1_bk[i], msd1_if[i]), file=anaout)


def getDomain(concentration, aP, aPcoords, Dxyz):
    if concentration == 'c000':
        domains = [29.387, 9.635, 97.136, 10.542]  # c000
    elif concentration == 'c020':
        domains = [28.330, 11.155, 100.320, 10.950]  # c020
    elif concentration == 'c040':
        domains = [26.760, 10.590, 102.400, 11.142]  # c040
    elif concentration == 'c060':
        domains = [26.140, 10.776, 105.340, 10.446]  # c060
    elif concentration == 'c080':
        domains = [25.320, 10.667, 107.740, 11.429]  # c080
    if1 = [domains[0] - domains[1] / 2, domains[0] + domains[1] / 2]
    bk1 = [domains[0] + domains[1] / 2, domains[2] - domains[3] / 2]
    if2 = [domains[2] - domains[3] / 2, domains[2] + domains[3] / 2]
    to1 = [domains[0] - domains[1] / 2, domains[2] + domains[3] / 2]
    aPinif1 = []
    aPinif2 = []
    aPinbk  = []
    aPinto  = []
    for i in range(len(aP)):
        if (aPcoords[:, Dxyz][i] >= if1[0] and aPcoords[:, Dxyz][i] <= if1[1]):
            aPinif1.append(i)
        if (aPcoords[:, Dxyz][i] >= if2[0] and aPcoords[:, Dxyz][i] <= if2[1]):
            aPinif2.append(i)
        if (aPcoords[:, Dxyz][i] >= bk1[0] and aPcoords[:, Dxyz][i] <= bk1[1]):
            aPinbk.append(i)
        if (aPcoords[:, Dxyz][i] >= to1[0] and aPcoords[:, Dxyz][i] <= to1[1]):
            aPinto.append(i)
    aPinif = aPinif1 + aPinif2
    return aPinto, aPinbk, aPinif


def main():
    # load trajectory for MSD
    concentration = 'c000'
    Dxyz = 2
    top = "../cg_topol.tpr"
    trj = "../submsd/unwrapped.xtc"
    trj2 = "../cg_traj.xtc"
    uta = mda.Universe(top, trj)
    uta2 = mda.Universe(top, trj2)
    aP = uta.select_atoms('type TF')  # fetch anions
    aP2 = uta2.select_atoms('type TF')  # fetch anions
    lP = np.linspace(0, len(aP) - 1, len(aP), dtype = int)  # create the atom index
    START = np.linspace(0, 8000, 201, dtype = int)
    tau = 2000
    MSDg0 = []
    MSDg1_to = []
    MSDg1_bk = []
    MSDg1_if = []
    for start in START:
        print(start)
        COORDs = []
        for ts in uta.trajectory:
            if ts.frame >= start:
                COORDs.append(aP.positions)  # unwrapped
            if ts.frame == start + tau + 1:
                break
        dispg0 = np.zeros((tau + 1, len(aP)))
        dispg1_to = np.zeros((tau + 1, len(aP)))
        dispg1_bk = np.zeros((tau + 1, len(aP)))
        dispg1_if = np.zeros((tau + 1, len(aP)))
        stateg1_to, stateg1_bk, stateg1_if = np.zeros(len(aP)), np.zeros(len(aP)), np.zeros(len(aP))
        firstframe = True
        frame = 0
        for time in range(start, tau + start + 1):
            if firstframe:  # marker the inital status of the ion group
                ts = uta2.trajectory[start]
                aPinto, aPinbk, aPinif = getDomain(concentration, aP2, aP2.positions, Dxyz)
                for atom in lP:
                    if atom in aPinbk:
                        stateg1_bk[atom] = 1
                    if atom in aPinif:
                        stateg1_if[atom] = 1
                    if atom in aPinto:
                        stateg1_to[atom] = 1
                firstframe = False
            displ = COORDs[frame] - COORDs[0]
            msd = np.sum(np.square(displ), axis = 1)
            dispg0[frame] = msd
            dispg1_to[frame] = msd * stateg1_to
            dispg1_bk[frame] = msd * stateg1_bk
            dispg1_if[frame] = msd * stateg1_if
            frame += 1
        msdg0 = np.zeros(tau + 1)
        msdg1_to = np.zeros(tau + 1)
        msdg1_bk = np.zeros(tau + 1)
        msdg1_if = np.zeros(tau + 1)
        for time in range(1, tau + 1):
            msdg0[time] = np.average(dispg0[time])
            msdg1_to[time] = dispg1_to[time][np.nonzero(dispg1_to[time])].mean()
            msdg1_bk[time] = dispg1_bk[time][np.nonzero(dispg1_bk[time])].mean()
            msdg1_if[time] = dispg1_if[time][np.nonzero(dispg1_if[time])].mean()
        MSDg0.append(msdg0)
        MSDg1_to.append(np.nan_to_num(msdg1_to))
        MSDg1_bk.append(np.nan_to_num(msdg1_bk))
        MSDg1_if.append(np.nan_to_num(msdg1_if))
    msd0 = np.zeros(tau + 1)
    msd1_to, msd1_bk, msd1_if = np.zeros(tau + 1), np.zeros(tau + 1), np.zeros(tau + 1)
    for time in range(1, tau + 1):
        tmp0 = []
        tmp1_to = []
        tmp1_bk = []
        tmp1_if = []
        for i in range(len(START)):
            tmp0.append(MSDg0[i][time])
            tmp1_to.append(MSDg1_to[i][time])
            tmp1_bk.append(MSDg1_bk[i][time])
            tmp1_if.append(MSDg1_if[i][time])
        tmp0 = np.asarray(tmp0)
        tmp1_to = np.asarray(tmp1_to)
        tmp1_bk = np.asarray(tmp1_bk)
        tmp1_if = np.asarray(tmp1_if)
        msd0[time] = np.average(tmp0)
        msd1_to[time] = tmp1_to[np.nonzero(tmp1_to)].mean()
        msd1_bk[time] = tmp1_bk[np.nonzero(tmp1_bk)].mean()
        msd1_if[time] = tmp1_if[np.nonzero(tmp1_if)].mean()
    xdata = np.linspace(0, tau, tau + 1)
    msdPrint('subMSD_domains', xdata, msd0, msd1_to, msd1_bk, msd1_if)

if __name__ == "__main__":
    main()