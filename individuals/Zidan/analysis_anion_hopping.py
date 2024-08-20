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

from analysis import *
import numpy as np
import seaborn as sn
import pickle
import networkx as nx
from collections import Counter


def hopPrint(filename, anion1, anion2, anion3, anion4):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# anion1, anion2, anion3, anion4', file=anaout)
        for i in range(len(anion1)):
            print('{:5d} {:5d} {:5d} {:5d}'.format(anion1[i], anion2[i], anion3[i], anion4[i]), file=anaout)


def hopes(assot, asso0):
    hopes = list(np.setdiff1d(assot, asso0))
    return hopes


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
    return set(at1), set(at2), type1, type2


def main():
    ac = h5py.File('ac.h5', 'r')
    al = h5py.File('al.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')
    END = 10000
    hop1, hop2, hop3, hop4 = [], [], [], []
    anion1, anion2, anion3, anion4 = [], [], [], []
    group1, group2 = [], []
    firstframe = True
    for time in range(0, END):
        t0 = time
        t1 = time + 1
        t2 = time + 2
        if t2 == END:
            break
        if firstframe:
            at1t0, at2t0, type1t0, type2t0 = getgroups(t0, ASSOac, ASSOal)
            at1t1, at2t1, type1t1, type2t1 = getgroups(t1, ASSOac, ASSOal)
            firstframe = False
        at1t2, at2t2, type1t2, type2t2 = getgroups(t2, ASSOac, ASSOal)
        hopes1 = hopes(type1t2, type1t1)
        hopes2 = hopes(type2t2, type2t1)
        hopt1, hopt2, hopts1, hopts2 = 0, 0, 0, 0
        at1t, at2t = [], []
        at1ts, at2ts = [], []
        for pair in hopes1:
            atom = pair[0]
            if atom in at1t1:
                at1t.append(atom)
            elif atom in at2t1 and pair not in type2t1:
                at1ts.append(atom)
            else:
                continue
        for pair in hopes2:
            atom = pair[0]
            if atom in at2t1:
                at2t.append(atom)
            elif atom in at1t1 and pair not in type1t1:
                at2ts.append(atom)
            else:
                continue
        anion1.append(len(set(at1t)))
        anion2.append(len(set(at2t)))
        anion3.append(len(set(at1ts)))
        anion4.append(len(set(at2ts)))
        at1t0, at2t0, type1t0, type2t0 = at1t1, at2t1, type1t1, type2t1
        at1t1, at2t1, type1t1, type2t1 = at1t2, at2t2, type1t2, type2t2
    hopPrint('group', anion1, anion2, anion3, anion4)


if __name__ == "__main__":
    main()