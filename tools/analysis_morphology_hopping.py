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
from datetime import datetime


def hopPrint(filename, anion1, anion2, anion3, anion4):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# anion1, anion2, anion3, anion4', file=anaout)
        for i in range(len(anion1)):
            print('{:5d} {:5d} {:5d} {:5d}'.format(anion1[i], anion2[i], anion3[i], anion4[i]), file=anaout)


def timestamp(TIME):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(TIME + " =", current_time)


def hopes(assot, asso0):
    hopes = list(np.setdiff1d(assot, asso0))
    return hopes


def getgroups(frame, ASSObk, ASSOif):
    abk = list(set([x[0] for x in ASSObk[frame]]))
    aif = list(set([x[0] for x in ASSOif[frame]]))
    return set(abk), set(aif), ASSObk[frame], ASSOif[frame]


def main():
    timestamp("Start Time")
    h5 = h5py.File('pils.h5', 'r')
    ASSObk = rawLoad(h5, 'in_htt')
    ASSOif = rawLoad(h5, 'it_htt')
    END = 10000
    anion1, anion2, anion3, anion4 = [], [], [], []
    firstframe = True
    for time in range(0, END):
        t0 = time
        t1 = time + 1
        if firstframe:
            at1t0, at2t0, type1t0, type2t0 = getgroups(t0, ASSObk, ASSOif)
            firstframe = False
        at1t1, at2t1, type1t1, type2t1 = getgroups(t1, ASSObk, ASSOif)
        hopes1 = hopes(type1t1, type1t0)
        hopes2 = hopes(type2t1, type2t0)
        at1t, at2t = [], []
        at1ts, at2ts = [], []
        for pair in hopes1:
            atom = pair[0]
            if atom in at1t0:
                at1t.append(atom)
            elif atom in at2t0 and pair not in type2t0:
                at1ts.append(atom)
            else:
                continue
        for pair in hopes2:
            atom = pair[0]
            if atom in at2t0:
                at2t.append(atom)
            elif atom in at1t0 and pair not in type1t0:
                at2ts.append(atom)
            else:
                continue
        anion1.append(len(set(at1t)))
        anion2.append(len(set(at2t)))
        anion3.append(len(set(at1ts)))
        anion4.append(len(set(at2ts)))
        at1t0, at2t0, type1t0, type2t0 = at1t1, at2t1, type1t1, type2t1
    hopPrint('group', anion1, anion2, anion3, anion4)
    timestamp("End Time")
    print("Hopping analysis finished!")


if __name__ == "__main__":
    main()