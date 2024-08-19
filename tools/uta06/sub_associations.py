"""
Copyright (C) 2018-2020 Zidan Zhang <zhangzidan@gmail.com>

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
import math
from collections import Counter

def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def assoFreq(frequence):
    """
    This function gets the normolized probability for association events.
    """
    fNum = frequence[0]
    for i in range(1, len(frequence)):
        fNum += frequence[i]
    numAsso = {}
    for k, v in fNum.items():
        numAsso[k] = v / float(len(frequence))
    totalAsso = sum(numAsso.values())
    avgAsso = {}
    for k, v in numAsso.items():
        avgAsso[k] = v / float(totalAsso)
    pairNum = []
    pairPro = []
    for k in sorted(avgAsso.keys()):
        pairNum.append(k)
        pairPro.append(avgAsso[k])
    return pairNum, pairPro


def anaPrint(filename, arrayNum, arrayPro):  # print association properties
    anaout = open(filename + '.dat', 'w')
    print('# ' + filename + ' Probability', file=anaout)
    for i in range(0, len(arrayNum)):
        print('{} {:10.5f}'.format(arrayNum[i], arrayPro[i]), file=anaout)
    anaout.close()


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


def getgroups(frame, ASSOac, ASSOal):
    poly = []; both = []
    ion = list(set([x[0] for x in ASSOal[frame]]))
    for i in range(len(ASSOac[frame])):
        polyIL = True  ## the given ion is only assiciated with the polymer
        atom = ASSOac[frame][i][0]
        if atom in ion:
            both.append(ASSOac[frame][i])
            polyIL = False
        if polyIL:
            poly.append(ASSOac[frame][i])
    return poly, both
        

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
    ac = h5py.File('../h5ac/pils.h5', 'r')
    al = h5py.File('../h5al/pils.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')
    freqAC_poly, freqCH_poly = [], []
    freqAC_both, freqCH_both = [], []
    frames = 10000
    for t0 in range(frames):
        poly, both = getgroups(t0, ASSOac, ASSOal)
        ac_poly, ch_poly = subasso(poly, 15)
        freqAC_poly.append(ac_poly)
        freqCH_poly.append(ch_poly)
        ac_both, ch_both = subasso(both, 15)
        freqAC_both.append(ac_both)
        freqCH_both.append(ch_both)
        if t0 % 100 == 0:
            print("Processing {:10.3f}".format((t0 + 1) * 1.0 / frames))
    acNum_poly, acPro_poly = assoFreq(freqAC_poly)
    chNum_poly, chPro_poly = assoFreq(freqCH_poly)
    anaPrint('polyIL_assoPair', acNum_poly, acPro_poly)
    anaPrint('polyIL_assoChain', chNum_poly, chPro_poly)
    acNum_both, acPro_both = assoFreq(freqAC_both)
    chNum_both, chPro_both = assoFreq(freqCH_both)
    anaPrint('Both_assoPair', acNum_both, acPro_both)
    anaPrint('Both_assoChain', chNum_both, chPro_both)


if __name__ == "__main__":
    main()
