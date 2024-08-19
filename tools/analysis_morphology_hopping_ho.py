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


def hopPrint(filename, anion):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# anion', file=anaout)
        for i in range(len(anion)):
            print('{:5d}'.format(anion[i]), file=anaout)


def timestamp(TIME):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(TIME + " =", current_time)


def hopes(assot, asso0):
    hopes = list(np.setdiff1d(assot, asso0))
    return hopes


def getgroups(frame, ASSO):
    anion = list(set([x[0] for x in ASSO[frame]]))
    return set(anion), ASSO[frame]


def main():
    timestamp("Start Time")
    h5 = h5py.File('pils.h5', 'r')
    ASSO = rawLoad(h5, 'htt')
    END = 10000
    anion1, anion2 = [], []
    firstframe = True
    for time in range(0, END):
        t0 = time
        t1 = time + 1
        if firstframe:
            at1t0, type1t0 = getgroups(t0, ASSO)
            firstframe = False
        at1t1, type1t1 = getgroups(t1, ASSO)
        hopes1 = hopes(type1t1, type1t0)
        at1t = [x[0] for x in hopes1]
        anion1.append(len(set(at1t)))
        at1t0, type1t0 = at1t1, type1t1
    hopPrint('group', anion1)
    timestamp("End Time")
    print("Hopping analysis finished!")


if __name__ == "__main__":
    main()