"""
Copyright (C) 2018-2021 Zidan Zhang <zhangzidan@gmail.com>

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


def groupPrint(filename, groups):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# len(at1), len(at2), len(at3), len(lt1), len(lt2), len(at31), len(at32)', file=anaout)
        for i in range(len(groups)):
            print(' '.join(str(x) for x in groups[i]), file=anaout)


def getgroups(frame, ASSOac, ASSOal):  # get the co-coordination anions
    at1, at2 = [], []
    anion = list(set([x[0] for x in ASSOal[frame]]))  # find the anions in co-coordination
    total = anion.copy()
    for pair in ASSOac[frame]:
        atom = pair[0]
        if atom not in anion:
            at1.append(atom)
        else:
            at2.append(atom)
        total.append(atom)
    at3 = set(total) - set(at1) - set(at2)
    return set(at1), set(at2), set(at3)


def ligroups(frame, ASSOac, ASSOla):
    at1, at2 = [], []
    coion = list(set([x[0] for x in ASSOac[frame]]))
    firstli = True
    for pair in ASSOla[frame]:
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
        if pair == ASSOla[frame][-1]:
            if set(liasso) & set(coion):
                at1.append(cur_li)
            else:
                at2.append(cur_li)
    return list(set(at1)), list(set(at2))


def subgroups(frame, ASSOac, ASSOal, at3, lig1, lig2):  # get the co-coordination anions
    at1, at2 = [], []
    firstion = True
    for pair in ASSOal[frame]:
        anion, lithium = pair
        if firstion:
            asso = []
            cur_ion = anion
            asso.append(lithium)
            firstion = False
        else:
            if anion == cur_ion:
                asso.append(lithium)
            else:
                if set(asso) & set(lig1) and anion in at3:
                    at1.append(cur_ion)
                elif set(asso) & set(lig2) and anion in at3:
                    at2.append(cur_ion)
                cur_ion = anion
                asso = []
                asso.append(lithium)
        if pair == ASSOal[frame][-1]:
            if set(asso) & set(lig1) and anion in at3:
                at1.append(cur_ion)
            elif set(asso) & set(lig2) and anion in at3:
                at2.append(cur_ion)
    return set(at1), set(at2)


def main():
    ac = h5py.File('../h5ac/pils.h5', 'r')
    al = h5py.File('../h5al/pils.h5', 'r')
    la = h5py.File('../h5la/pils.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')
    ASSOla = rawLoad(la, 'htt')
    groups = []
    frames = 10001
    for time in np.linspace(0, frames - 1, 1001, dtype = int):
        if time % 100 == 0:
            print("Processing {:10.3f}".format(time))
        at1, at2, at3 = getgroups(time, ASSOac, ASSOal)
        lt1, lt2 = ligroups(time, ASSOac, ASSOla)
        at31, at32 = subgroups(time, ASSOac, ASSOal, at3, lt1, lt2)
        groups.append([len(at1), len(at2), len(at3), len(lt1), len(lt2), len(at31), len(at32)])
    groupPrint('ionFraction', groups)


if __name__ == "__main__":
    main()