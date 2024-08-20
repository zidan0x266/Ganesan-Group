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
from collections import Counter


def hopPrint(filename, anion1, anion2, anion3, anion21, anion31, anion12, anion32, anion13, anion23):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# anion1, anion2, anion3, anion21, anion31, anion12, anion32, anion13, anion23', file=anaout)
        for i in range(len(anion1)):
            print('{:5d} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d}'.format(anion1[i], anion2[i], anion3[i], anion21[i], anion31[i], anion12[i], anion32[i], anion13[i], anion23[i]), file=anaout)


def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def hopes(assot, asso0):
    hopes = list(np.setdiff1d(assot, asso0))
    return hopes


def getgroups(frame, lP, ASSOac, ASSOal):  # get the co-coordination anions
    type1, type2 = [], []
    type3, type4 = [], []
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
    atomg1 = set(at1)
    atomg2 = set(at2)
    atomg3 = np.setdiff1d(lP, list(atomg1) + list(atomg2))
    for pair in ASSOal[frame]:
        atom = pair[0]
        if atom in atomg3:
            type3.append(pair)
        else:
            type4.append(pair)
    return atomg1, atomg2, atomg3, type1, type2, type3, type4


def main():
    top = "../md.tpr"
    uta = mda.Universe(top) 
    aP = uta.select_atoms('type opls_480')  # fetch anions
    lP = np.linspace(0, len(aP) - 1, len(aP), dtype = int)  # create the atom index
    ac = h5py.File('../h5ac/pils.h5', 'r')
    al = h5py.File('../h5al/pils.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')
    END = 10001
    anion1, anion2, anion3 = [], [], []
    anion21, anion31 = [], []
    anion12, anion32 = [], []
    anion13, anion23 = [], []
    STEP = 1
    for time in range(0, END):
        if time % 100 == 0:
            print("Processing {:10.3f}".format((time + 1) * 1.0 / END))
        t0 = time
        t1 = time + STEP
        if t1 == END - 1:
            break
        at1t0, at2t0, at3t0, type1t0, type2t0, type3t0, type4t0 = getgroups(t0, lP, ASSOac, ASSOal)
        at1t1, at2t1, at3t1, type1t1, type2t1, type3t1, type4t1 = getgroups(t1, lP, ASSOac, ASSOal)
        hopes1 = hopes(type1t1, type1t0)
        hopes2 = hopes(type2t1, type2t0)
        hopes3 = hopes(type3t1, type3t0)
        hopes4 = hopes(type4t1, type4t0)
        at1t, at2t, at3t = [], [], []
        ats21, ats31 = [], []
        ats12, ats32 = [], []
        ats13, ats23 = [], []
        for pair in hopes1:
            atom = pair[0]
            if atom in at1t0:
                at1t.append(atom)
            elif atom in at2t0:
                ats21.append(atom)
            elif atom in at3t0:
                ats31.append(atom)
            else:
                continue
        for pair in hopes2:
            atom = pair[0]
            if atom in at2t0:
                at2t.append(atom)
            elif atom in at1t0:
                ats12.append(atom)
            elif atom in at3t0:
                ats32.append(atom)
            else:
                continue
        for pair in hopes3:
            atom = pair[0]
            if atom in at3t0:
                at3t.append(atom)
            elif atom in at1t0:
                ats13.append(atom)
            elif atom in at2t0:
                ats23.append(atom)
            else:
                continue
        for pair in hopes4:
            atom = pair[0]
            if atom in at2t0:
                at2t.append(atom)
            elif atom in at1t0:
                ats12.append(atom)
            elif atom in at3t0:
                ats32.append(atom)
            else:
                continue
        anion1.append(len(set(at1t)))
        anion2.append(len(set(at2t)))
        anion3.append(len(set(at3t)))
        anion21.append(len(set(ats21)))
        anion31.append(len(set(ats31)))
        anion12.append(len(set(ats12)))
        anion32.append(len(set(ats32)))
        anion13.append(len(set(ats13)))
        anion23.append(len(set(ats23)))
    hopPrint("switching_$STEP", anion1, anion2, anion3, anion21, anion31, anion12, anion32, anion13, anion23)


if __name__ == "__main__":
    main()