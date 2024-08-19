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
from analysis import *

def main():
    ac = h5py.File('ac.h5', 'r')
    al = h5py.File('al.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')
    anion1, anion2, anion3 = [], [], []
    for frame in range(len(ASSOac)):
        at1, at2 = [], []
        cation = []
        lithium = []
        for i in range(len(ASSOac[frame])):
            atom1, atom2 = ASSOac[frame][i]
            at1.append(atom1)
            cation.append(atom2)
        for i in range(len(ASSOal[frame])):
            atom1, atom2 = ASSOal[frame][i]
            at2.append(atom1)
            lithium.append(atom2)
        anion1.append(len(list(np.setdiff1d(at1, at2))))
        anion2.append(len(list(np.setdiff1d(at2, at1))))
        anion3.append(360 - len(list(np.setdiff1d(at1, at2))) - len(list(np.setdiff1d(at2, at1))))
    print(np.average(anion1), np.average(anion1), np.average(anion1))

if __name__ == "__main__":
    main()