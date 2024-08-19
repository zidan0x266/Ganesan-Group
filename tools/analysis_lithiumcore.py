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
from collections import Counter
from analysis import *

def outPrint(filename, xdata, ydata, zdata):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# jump1 jump2 jump3', file=anaout)
        for i in range(len(xdata)):
            print('{:5d} {:5d} {:5d}'.format(xdata[i], ydata[i], zdata[i]), file=anaout)
        

def LithiumCore(ASSO, LN):
    LCore = []
    for lithium in range(LN):
        for i in range(len(ASSO)):
            atomi, atomj = ASSO[i]
            if atomj == lithium:
                LCore.append((atomj, atomi))
    return LCore

    

def subasso(ASSO):
    """
    This function outputs the analysis for association relationship of 
    a given sub-domain.

    Process a single frame.
    """
    assoAC = []  # total associations
    firstion = True
    for i in range(len(ASSO)):
        if firstion:
            t1 = 0
            atom0 = ASSO[0][0]
            firstion = False
        atom1, atom2 = ASSO[i]
        if atom1 == atom0:
            t1 += 1
        else:
            assoAC.append(t1)
            t1 = 0
            t1 += 1
            atom0 = atom1
    return Counter(assoAC)


def main():
    al = h5py.File('pils.h5', 'r')
    ASSOal = rawLoad(al, 'htt')
    freqAC = []
    for t in range(100):
        ASSOla = LithiumCore(ASSOal[t], 30)
        ac = subasso(ASSOla)
        freqAC.append(ac)
        if (t + 1) % 10 == 0:
            print("Processing {:10.3f}".format((t + 1) / 100))
    acNum, acPro = assoFreq(freqAC)
    anaPrint('LithiumCore_assoPair', acNum, acPro)


if __name__ == "__main__":
    main()