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

def outPrint(filename, xdata, ydata, zdata):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# jump1 jump2 jump3', file=anaout)
        for i in range(len(xdata)):
            print('{:5d} {:5d} {:5d}'.format(xdata[i], ydata[i], zdata[i]), file=anaout)


def intintPrint(filename, xdata, ydata):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# Intermolecular Intramolecular', file=anaout)
        for i in range(len(xdata)):
            print('{:5d} {:5d}'.format(xdata[i], ydata[i]), file=anaout)


def getchain(cocoset, DP):
    firstli = True
    linfo = []
    lichain = []
    t1 = 0
    for coco in cocoset:
        if firstli:
            prevli = coco[0]
            firstli = False
        lidx = coco[0]  # lithium index
        if lidx == prevli:
            t1 += 1
            catdx = coco[2]  # polycation index
            chaindx = int(math.floor(catdx / DP) + 1)  # polymer chain index
            linfo.append((lidx, chaindx))
        else:
            lichain.append(t1)
            t1 = 1
            prevli = lidx
            catdx = coco[2]  # polycation index
            chaindx = int(math.floor(catdx / DP) + 1)  # polymer chain index
            linfo.append((lidx, chaindx))
    linform = sorted(list(set(linfo)), key=lambda lion: lion[0])
    return Counter(lichain), linform

def getcoco(frame, ASSOac, ASSOal):  # the cation is the reference
    coco = []
    ion =[]
    for pair in range(len(ASSOal[frame])):
        ion.append(ASSOal[frame][pair][0])
    for i in range(len(ASSOac[frame])):
        atom1, atom2 = ASSOac[frame][i]
        if atom1 in ion:
            for j in range(len(ASSOal[frame])):
                atom3, atom4 = ASSOal[frame][j]
                if atom3 == atom1:
                    coco.append((atom2, atom1, atom4))
    return coco


def licoco(frame, ASSOac, ASSOal):  # the lithium is the reference
    coco = []
    ion =[]
    for pair in range(len(ASSOal[frame])):
        ion.append(ASSOal[frame][pair][0])
    for i in range(len(ASSOac[frame])):
        atom1, atom2 = ASSOac[frame][i]
        if atom1 in ion:
            for j in range(len(ASSOal[frame])):
                atom3, atom4 = ASSOal[frame][j]
                if atom3 == atom1:
                    coco.append((atom4, atom1, atom2))
    return coco


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


def main():
    ac = h5py.File('ac.h5', 'r')
    al = h5py.File('al.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')
    DP = 15
    chainfreq = []
    firstframe = True
    dinter, dintra = [], []
    for frame in range(0, 50):
        t = frame + 1
        if firstframe:
            cocot0 = licoco(frame, ASSOac, ASSOal)
            scoco0 = sorted(cocot0, key=lambda lion: lion[0])
            firstframe = False
        cocot1 = licoco(t, ASSOac, ASSOal)
        scoco1 = sorted(cocot1, key=lambda lion: lion[0])
        asso = []
        for coco in scoco1:
            if coco not in scoco0:
                asso.append(coco)
        chaindis, linfo = getchain(scoco0, DP)
        chainfreq.append(chaindis)
        chaindis1, linfo1 = getchain(asso, DP)
        scoco0 = scoco1
        inter, intra = 0, 0
        for x in linfo1:
            if x in linfo:
                inter += 1
            else:
                intra += 1
        dinter.append(inter)
        dintra.append(intra)
    acNum, acPro = assoFreq(chainfreq)
    anaPrint('LiChain', acNum, acPro)
    intintPrint('LiInterIntra', dinter, dintra)


if __name__ == "__main__":
    main()
