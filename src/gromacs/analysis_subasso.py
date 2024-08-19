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


def getcoco(frame, ASSOac, ASSOal):
    coco = []
    for i in range(len(ASSOac[frame])):
        atom1, atom2 = ASSOac[frame][i]
        for j in range(len(ASSOal[frame])):
            atom3, atom4 = ASSOal[frame][j]
            if atom3 == atom1:
                coco.append((atom2, atom1, atom4))
    return coco


def calcNode(ASSOac):
    asso_connections = set()
    asso_activated = set()
    for i in range(len(ASSOac)):
        atom1, atom2 = ASSOac[i]
        asso_connections.add(tuple(sorted((atom1, atom2))))
        asso_activated.add(atom1)
        asso_activated.add(atom2)
    return asso_activated, asso_connections


def getpoly(frame, ASSOac, ASSOal):
    poly = []
    ion = []
    for pair in range(len(ASSOal[frame])):
        ion.append(ASSOal[frame][pair][0])
    for i in range(len(ASSOac[frame])):
        polyIL = True  ## the given ion is only assiciated with the polymer
        atom = ASSOac[frame][i][0]
        if atom in ion:
            polyIL = False
        if polyIL:
            poly.append(ASSOac[frame][i])
    return poly


def getboth(frame, ASSOac, ASSOal):
    both = []
    ion = []
    for pair in range(len(ASSOal[frame])):
        ion.append(ASSOal[frame][pair][0])
    for i in range(len(ASSOac[frame])):
        atom = ASSOac[frame][i][0]
        if atom in ion:
            both.append(ASSOac[frame][i])
    return both


def hopes(coco1, coco0):
    hopes, polyIL, lithium = [], [], []
    counter = 0
    for i in coco1:
        if i not in coco0:
            counter += 1
            hopes.append(i)
    return hopes, counter


def cococounter(hopes, polyIL0, lithium0):
    jump1, jump2, jump3 = 0, 0, 0
    base1, base2 = [], []
    for i in range(len(polyIL0)):
        base1.append((polyIL0[i][0], polyIL0[i][1]))
    for i in range(len(lithium0)):
        base2.append((lithium0[i][0], lithium0[i][1]))
    for i in hopes:
        at1, at2, at3 = i[0], i[1], i[2]
        polyIL = (at2, at1)
        lithium = (at2, at3)
        if polyIL not in base1 and lithium not in base2:
            jump3 += 1
        elif polyIL in base1 and lithium not in base2:
            jump1 += 1
        elif polyIL not in base1 and lithium in base2:
            jump2 += 1
    return jump1, jump2, jump3
        

def subasso(ASSO, counterN):
    """
    This function outputs the analysis for association relationship of 
    a given sub-domain.

    Process a single frame.
    """
    assoAC = []; assoCH = [] # total associations
    firstion = True
    for i in range(len(ASSO)):
        if firstion:
            beads = []; chains = []
            t1 = 0
            atom0 = ASSO[0][0]
            firstion = False
        atom1, atom2 = ASSO[i]
        if atom1 == atom0:
            t1 += 1
            beads.append(atom2)
        else:
            if len(beads) != 0:
                for k in range(0, len(beads)):
                    t3 = int(math.floor(beads[k] / counterN) + 1)
                    chains.append(t3)
            else:
                chains.append(0)
            if len(chains) == 1 and chains[0] == 0:
                assoCH.append(0)
            else:
                assoCH.append(len(Counter(chains)))
            assoAC.append(t1)
            beads = []; chains = []
            t1 = 0
            t1 += 1
            beads.append(atom2)
            atom0 = atom1
    return Counter(assoAC), Counter(assoCH)

def main():
    ac = h5py.File('ac.h5', 'r')
    al = h5py.File('al.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')
    freqAC, freqCH = [], []
    for t0 in range(100):
        poly = getpoly(t0, ASSOac, ASSOal)
        ac, ch = subasso(poly, 15)
        freqAC.append(ac)
        freqCH.append(ch)
        if t0 % 100 == 0:
            print("Processing {:10.3f}".format(t0 / 2000))
    acNum, acPro = assoFreq(freqAC)
    chNum, chPro = assoFreq(freqCH)
    anaPrint('polyIL_assoPair', acNum, acPro)
    anaPrint('polyIL_assoChain', chNum, chPro)
    freqAC, freqCH = [], []
    for t0 in range(100):
        both = getboth(t0, ASSOac, ASSOal)
        ac, ch = subasso(both, 15)
        freqAC.append(ac)
        freqCH.append(ch)
        if t0 % 100 == 0:
            print("Processing {:10.3f}".format(t0 / 2000))
    acNum, acPro = assoFreq(freqAC)
    chNum, chPro = assoFreq(freqCH)
    anaPrint('Both_assoPair', acNum, acPro)
    anaPrint('Both_assoChain', chNum, chPro)


if __name__ == "__main__":
    main()