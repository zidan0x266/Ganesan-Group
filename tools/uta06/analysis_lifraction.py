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


def ligroups(frame, ASSOac, ASSOal):
    at1, at2 = [], []
    coion = list(set([x[0] for x in ASSOac[frame]]))
    firstli = True
    for pair in ASSOal[frame]:
        lithium, anion = pair
        if firstli:
            liasso = []
            cur_li = lithium
            liasso.append(anion)
            firstli = False
        else:
            if lithium == cur_li:
                liasso.append(anion)
            else:
                if set(liasso) & set(coion):
                    at1.append(cur_li)
                else:
                    at2.append(cur_li)
                cur_li = lithium
                liasso = []
                liasso.append(anion)
        if pair == ASSOal[frame][-1]:
            if set(liasso) & set(coion):
                at1.append(cur_li)
            else:
                at2.append(cur_li)
    return list(set(at1)), list(set(at2))


def main():
    ac = h5py.File('../h5ac/pils.h5', 'r')
    al = h5py.File('../h5la/pils.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')
    polyIL, Only = [], []
    frames = 10001
    for time in range(frames):
        if time % 100 == 0:
            print("Processing {:10.3f}".format((time + 1) * 1.0 / frames))
        poly, only = ligroups(time, ASSOac, ASSOal)
        polyIL.append(len(poly))
        Only.append(len(only))
    groupPrint('lithiumFraction', polyIL, Only)


if __name__ == "__main__":
    main()
