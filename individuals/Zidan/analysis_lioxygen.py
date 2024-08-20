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
from collections import Counter

def getoxygen(ASSO, DP):
    firstli = True
    litfsi = []
    t1 = 0
    for asso in ASSO:
        lidx = asso[0]  # lithium index
        oxydx = asso[1]  # polycation index
        tfsi = int(math.floor(oxydx / DP) + 1)  # polymer chain index
        litfsi.append((lidx, tfsi))
    mono, bi = 0, 0
    for _, val in Counter(litfsi).items():
        if val == 1:
            mono += 1
        if val == 2:
            bi += 1
    return mono, bi


def main():
    conc = 'c040'
    samples = ['s01', 's02', 's03', 's04', 's05', 
              's06', 's07', 's08', 's09', 's10']
    mono, bi = [], []
    for sample in samples:
        lo = h5py.File('../' + sample + '/' + conc + '/' + 'analysis_lo/pils.h5', 'r')
        for frame in range(10000):
            ASSOlo = rawLoad(lo, 'htt')[frame]
            x, y = getoxygen(ASSOlo, 4)
            mono.append(x)
            bi.append(y)
    print(np.average(mono), np.average(bi))
    print(np.std(mono), np.std(bi))


if __name__ == "__main__":
    main()
