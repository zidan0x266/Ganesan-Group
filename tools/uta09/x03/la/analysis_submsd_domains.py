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


def msdPrint(filename, time, msd0, msd1_to, msd2_to, msd3_to, msd1_bk, msd2_bk, msd3_bk, msd1_if, msd2_if, msd3_if):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# time, total, msd1_to, msd2_to, msd3_to, msd1_bk, msd2_bk, msd3_bk, msd1_if, msd2_if, msd3_if', file=anaout)
        for i in range(len(time)):
            print('{:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(time[i], msd0[i], 
                                                                        msd1_to[i], msd2_to[i], msd3_to[i], 
                                                                        msd1_bk[i], msd2_bk[i], msd3_bk[i], 
                                                                        msd1_if[i], msd2_if[i], msd3_if[i]), file=anaout)


def getgroups(frame, ASSOac, ASSOal):  # get the co-coordination anions
    at1, at2 = [], []
    anion = list(set([x[0] for x in ASSOal[frame]]))  # find the anions in co-coordination
    for pair in ASSOac[frame]:
        atom = pair[0]
        if atom not in anion:
            at1.append(atom)
        else:
            at2.append(atom)
    return list(set(at1)), list(set(at2))


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
    concentration = 'c020'
    Dxyz = 0
    top = "../cg_topol.tpr"
    trj = "../travel/unwrapped.xtc"
    trj2 = "../cg_traj.xtc"
    uta = mda.Universe(top, trj)
    uta2 = mda.Universe(top, trj2)
    aP = uta.select_atoms('type TF')  # fetch anions
    aP2 = uta2.select_atoms('type TF')  # fetch anions
    lP = np.linspace(0, len(aP) - 1, len(aP), dtype = int)  # create the atom index
    START = np.linspace(0, 8000, 201, dtype = int)
    tau = 2000
    # load ion pair associations
    ac = h5py.File('../h5ac/pils.h5', 'r')
    al = h5py.File('../h5al/pils.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')
    MSDg0 = []
    MSDg1_to, MSDg2_to, MSDg3_to = [], [], []
    MSDg1_bk, MSDg2_bk, MSDg3_bk = [], [], []
    MSDg1_if, MSDg2_if, MSDg3_if = [], [], []
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
        dispg2_to = np.zeros((tau + 1, len(aP)))
        dispg3_to = np.zeros((tau + 1, len(aP)))
        dispg1_bk = np.zeros((tau + 1, len(aP)))
        dispg2_bk = np.zeros((tau + 1, len(aP)))
        dispg3_bk = np.zeros((tau + 1, len(aP)))
        dispg1_if = np.zeros((tau + 1, len(aP)))
        dispg2_if = np.zeros((tau + 1, len(aP)))
        dispg3_if = np.zeros((tau + 1, len(aP)))
        stateg1_to, stateg1_bk, stateg1_if = np.zeros(len(aP)), np.zeros(len(aP)), np.zeros(len(aP))
        stateg2_to, stateg2_bk, stateg2_if = np.zeros(len(aP)), np.zeros(len(aP)), np.zeros(len(aP))
        stateg3_to, stateg3_bk, stateg3_if = np.zeros(len(aP)), np.zeros(len(aP)), np.zeros(len(aP))
        firstframe = True
        frame = 0
        for time in range(start, tau + start + 1):
            if firstframe:  # marker the inital status of the ion group
                ts = uta2.trajectory[start]
                aPinto, aPinbk, aPinif = getDomain(concentration, aP2, aP2.positions, Dxyz)
                atomg1, atomg2 = getgroups(time, ASSOac, ASSOal)
                atomg3 = np.setdiff1d(lP, atomg1 + atomg2)  # ignore this contribution
                for atom in atomg1:
                    if atom in aPinbk:
                        stateg1_bk[atom] = 1
                    if atom in aPinif:
                        stateg1_if[atom] = 1
                    if atom in aPinto:
                        stateg1_to[atom] = 1
                for atom in atomg2:
                    if atom in aPinbk:
                        stateg2_bk[atom] = 1
                    if atom in aPinif:
                        stateg2_if[atom] = 1
                    if atom in aPinto:
                        stateg2_to[atom] = 1
                for atom in atomg3:
                    if atom in aPinbk:
                        stateg3_bk[atom] = 1
                    if atom in aPinif:
                        stateg3_if[atom] = 1
                    if atom in aPinto:
                        stateg3_to[atom] = 1
                firstframe = False
            displ = COORDs[frame] - COORDs[0]
            msd = np.sum(np.square(displ), axis = 1)
            dispg0[frame] = msd
            dispg1_to[frame] = msd * stateg1_to
            dispg2_to[frame] = msd * stateg2_to
            dispg3_to[frame] = msd * stateg3_to
            dispg1_bk[frame] = msd * stateg1_bk
            dispg2_bk[frame] = msd * stateg2_bk
            dispg3_bk[frame] = msd * stateg3_bk
            dispg1_if[frame] = msd * stateg1_if
            dispg2_if[frame] = msd * stateg2_if
            dispg3_if[frame] = msd * stateg3_if
            frame += 1
        msdg0 = np.zeros(tau + 1)
        msdg1_to = np.zeros(tau + 1)
        msdg2_to = np.zeros(tau + 1)
        msdg3_to = np.zeros(tau + 1)
        msdg1_bk = np.zeros(tau + 1)
        msdg2_bk = np.zeros(tau + 1)
        msdg3_bk = np.zeros(tau + 1)
        msdg1_if = np.zeros(tau + 1)
        msdg2_if = np.zeros(tau + 1)
        msdg3_if = np.zeros(tau + 1)
        for time in range(1, tau + 1):
            msdg0[time] = np.average(dispg0[time])
            msdg1_to[time] = dispg1_to[time][np.nonzero(dispg1_to[time])].mean()
            msdg2_to[time] = dispg2_to[time][np.nonzero(dispg2_to[time])].mean()
            msdg3_to[time] = dispg3_to[time][np.nonzero(dispg3_to[time])].mean()
            msdg1_bk[time] = dispg1_bk[time][np.nonzero(dispg1_bk[time])].mean()
            msdg2_bk[time] = dispg2_bk[time][np.nonzero(dispg2_bk[time])].mean()
            msdg3_bk[time] = dispg3_bk[time][np.nonzero(dispg3_bk[time])].mean()
            msdg1_if[time] = dispg1_if[time][np.nonzero(dispg1_if[time])].mean()
            msdg2_if[time] = dispg2_if[time][np.nonzero(dispg2_if[time])].mean()
            msdg3_if[time] = dispg3_if[time][np.nonzero(dispg3_if[time])].mean()
        MSDg0.append(msdg0)
        MSDg1_to.append(np.nan_to_num(msdg1_to))
        MSDg2_to.append(np.nan_to_num(msdg2_to))
        MSDg3_to.append(np.nan_to_num(msdg3_to))
        MSDg1_bk.append(np.nan_to_num(msdg1_bk))
        MSDg2_bk.append(np.nan_to_num(msdg2_bk))
        MSDg3_bk.append(np.nan_to_num(msdg3_bk))
        MSDg1_if.append(np.nan_to_num(msdg1_if))
        MSDg2_if.append(np.nan_to_num(msdg2_if))
        MSDg3_if.append(np.nan_to_num(msdg3_if))
    msd0 = np.zeros(tau + 1)
    msd1_to, msd1_bk, msd1_if = np.zeros(tau + 1), np.zeros(tau + 1), np.zeros(tau + 1)
    msd2_to, msd2_bk, msd2_if = np.zeros(tau + 1), np.zeros(tau + 1), np.zeros(tau + 1)
    msd3_to, msd3_bk, msd3_if = np.zeros(tau + 1), np.zeros(tau + 1), np.zeros(tau + 1)
    for time in range(1, tau + 1):
        tmp0 = []
        tmp1_to, tmp2_to, tmp3_to = [], [], []
        tmp1_bk, tmp2_bk, tmp3_bk = [], [], []
        tmp1_if, tmp2_if, tmp3_if = [], [], []
        for i in range(len(START)):
            tmp0.append(MSDg0[i][time])
            tmp1_to.append(MSDg1_to[i][time])
            tmp2_to.append(MSDg2_to[i][time])
            tmp3_to.append(MSDg3_to[i][time])
            tmp1_bk.append(MSDg1_bk[i][time])
            tmp2_bk.append(MSDg2_bk[i][time])
            tmp3_bk.append(MSDg3_bk[i][time])
            tmp1_if.append(MSDg1_if[i][time])
            tmp2_if.append(MSDg2_if[i][time])
            tmp3_if.append(MSDg3_if[i][time])
        tmp0 = np.asarray(tmp0)
        tmp1_to = np.asarray(tmp1_to)
        tmp2_to = np.asarray(tmp2_to)
        tmp3_to = np.asarray(tmp3_to)
        tmp1_bk = np.asarray(tmp1_bk)
        tmp2_bk = np.asarray(tmp2_bk)
        tmp3_bk = np.asarray(tmp3_bk)
        tmp1_if = np.asarray(tmp1_if)
        tmp2_if = np.asarray(tmp2_if)
        tmp3_if = np.asarray(tmp3_if)
        msd0[time] = np.average(tmp0)
        msd1_to[time] = tmp1_to[np.nonzero(tmp1_to)].mean()
        msd2_to[time] = tmp2_to[np.nonzero(tmp2_to)].mean()
        msd3_to[time] = tmp3_to[np.nonzero(tmp3_to)].mean()
        msd1_bk[time] = tmp1_bk[np.nonzero(tmp1_bk)].mean()
        msd2_bk[time] = tmp2_bk[np.nonzero(tmp2_bk)].mean()
        msd3_bk[time] = tmp3_bk[np.nonzero(tmp3_bk)].mean()
        msd1_if[time] = tmp1_if[np.nonzero(tmp1_if)].mean()
        msd2_if[time] = tmp2_if[np.nonzero(tmp2_if)].mean()
        msd3_if[time] = tmp3_if[np.nonzero(tmp3_if)].mean()
    xdata = np.linspace(0, tau, tau + 1)
    msdPrint('subMSD_domains', xdata, msd0, msd1_to, msd2_to, msd3_to, msd1_bk, msd2_bk, msd3_bk, msd1_if, msd2_if, msd3_if)

if __name__ == "__main__":
    main()