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


def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def groupPrint(filename, anion1, anion2):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# polyIL Both', file=anaout)
        for i in range(len(anion1)):
            print('{:5d} {:5d}'.format(anion1[i], anion2[i]), file=anaout)


def cationgroups(frame, ASSOac, ASSOal):
    at1, at2 = [], []
    coion = list(set([x[0] for x in ASSOal[frame]]))
    firstcation = True
    for pair in ASSOac[frame]:
        cation, anion = pair
        if firstcation:
            cationasso = []
            cur_cation = cation
            cationasso.append(anion)
            firstcation = False
        else:
            if cation == cur_cation:
                cationasso.append(anion)
            else:
                if set(cationasso) & set(coion):
                    at1.append(cur_cation)
                else:
                    at2.append(cur_cation)
                cur_cation = cation
                cationasso = []
                cationasso.append(anion)
        if pair == ASSOac[frame][-1]:
            if set(cationasso) & set(coion):
                at1.append(cur_cation)
            else:
                at2.append(cur_cation)
    return list(set(at1)), list(set(at2))


def main():
    ac = h5py.File('../h5ca/pils.h5', 'r')
    al = h5py.File('../h5al/pils.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')
    polyIL, Only = [], []
    frames = 10001
    for time in range(frames):
        if time % 100 == 0:
            print("Processing {:10.3f}".format((time + 1) * 1.0 / frames))
        poly, only = cationgroups(time, ASSOac, ASSOal)
        polyIL.append(len(poly))
        Only.append(len(only))
    groupPrint('cationFraction', polyIL, Only)


if __name__ == "__main__":
    main()