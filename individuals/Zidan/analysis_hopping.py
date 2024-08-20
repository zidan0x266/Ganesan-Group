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

import math
import numpy as np
import h5py
from collections import Counter

def outPrint(filename, xdata, ydata):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# vehicular structural', file=anaout)
        for i in range(len(xdata)):
            print('{:5d} {:5d}'.format(xdata[i], ydata[i]), file=anaout)


def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def getcoco(frame, ASSOac, ASSOal):
    coco = []
    anion = list(set([x[0] for x in ASSOal[frame]]))
    base1 = [(pair[0], pair[1]) for pair in ASSOac[frame]]
    base2 = [(pair[0], pair[1]) for pair in ASSOal[frame]]
    for i in base1:
        atom1, atom2 = i
        if atom1 in anion:
            for j in base2:
                atom3, atom4 = j
                if atom3 == atom1:
                    coco.append((atom2, atom1, atom4))
    return coco


def hopes(coco1, coco0):
    hopes = [item for item in coco1 if item not in coco0]
    return hopes


def cococounter(hopes, polyIL0, lithium0):
    base1 = [(pair[0], pair[1]) for pair in polyIL0]
    base2 = [(pair[0], pair[1]) for pair in lithium0]
    pair0, pair1, pair2 = [], [], []
    for i in hopes:
        at1, at2, at3 = i[0], i[1], i[2]
        polyIL = (at2, at1)
        lithium = (at2, at3)
        if polyIL not in base1 and lithium not in base2:
            pair0.append(lithium[1])
        elif polyIL in base1 and lithium not in base2:
            pair1.append(lithium[1])
        elif polyIL not in base1 and lithium in base2:
            pair2.append(lithium[1])
    jump1 = len(set(pair2) - set(pair1) - set(pair0))
    jump2 = len(set(pair0 + pair1 + pair2)) - jump1
    return jump1, jump2
        

def main():
    ac = h5py.File('ac.h5', 'r')
    al = h5py.File('al.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')
    jump1, jump2 = [], []
    firstframe = True
    for t0 in range(10000):
        t = t0 + 1
        if firstframe:
            cocot0 = getcoco(t0, ASSOac, ASSOal)
            firstframe = False
        cocot = getcoco(t, ASSOac, ASSOal)
        hops = hopes(cocot, cocot0)
        x, y = cococounter(hops, ASSOac[t0], ASSOal[t0])
        jump1.append(x)
        jump2.append(y)
        cocot0 = cocot
        print(x, y)
    outPrint("cocohop", jump1, jump2)


if __name__ == "__main__":
    main()