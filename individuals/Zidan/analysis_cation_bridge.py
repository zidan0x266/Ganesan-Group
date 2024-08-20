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
import math
from itertools import combinations

def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def anaPrint(filename, d1, d2, d3, d4):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# Btg1, Bsg1, Btg2, Bsg2', file=anaout)
        for i in range(len(d1)):
            print('{:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(d1[i], d2[i], d3[i], d4[i]), file=anaout)


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


def rSubset(arr, r):
    # return list of all subsets of length r 
    # to deal with duplicate subsets use  
    # set(list(combinations(arr, r))) 
    return list(combinations(arr, r))


def getbridges(ASSO):
    firstion = True
    bridges = []
    for i in range(len(ASSO)):
        if firstion:
            beads = []
            atom0 = ASSO[0][0]
            firstion = False
        atom1, atom2 = ASSO[i]
        if atom1 == atom0:
            beads.append(atom2)
        else:
            if len(beads) > 1:
                catcat = rSubset(beads, 2)
                bridges.append(catcat)
            beads = []
            beads.append(atom2)
            atom0 = atom1
    return bridges


def getarchs(BRIDGES, DP):
    bsingle = []
    btwins = []
    for i in range(len(BRIDGES)):
        bs, bt = 0, 0
        for pair in BRIDGES[i]:
            cat1, cat2 = pair
            t1 = int(math.floor(cat1 / DP) + 1)
            t2 = int(math.floor(cat2 / DP) + 1)
            if t1 != t2:
                bt += 1
            else:
                bs += 1
        bsingle.append(bs)
        btwins.append(bt)
    return bsingle, btwins


def main():
    # load ion pair associations
    ac = h5py.File('ac.h5', 'r')
    al = h5py.File('al.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')   
    Bsingle1, Btwins1 = [], []
    Bsingle2, Btwins2 = [], []
    for i in range(101):
        if i % 100 == 0:
            print("Currently processing the {} frame...".format(i))
        n1, n2, g1, g2 = getgroups(i, ASSOac, ASSOal)
        b1 = getbridges(g1)
        b2 = getbridges(g2)
        bs1, bt1 = getarchs(b1, 15)
        bs2, bt2 = getarchs(b2, 15)
        Btwins1.append(np.sum(bt1) / len(n1))
        Bsingle1.append(np.sum(bs1) / len(n1))
        Btwins2.append(np.sum(bt2) / len(n2))
        Bsingle2.append(np.sum(bs2) / len(n2))
    anaPrint('Bridge', Btwins1, Bsingle1, Btwins2, Bsingle2)
    print(np.average(Btwins1), np.std(Btwins1))
    print(np.average(Btwins2), np.std(Btwins2))
    print(np.average(Bsingle1), np.std(Bsingle1))
    print(np.average(Bsingle2), np.std(Bsingle2))

if __name__ == "__main__":
    main()