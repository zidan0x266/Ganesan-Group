"""
Copyright (C) 2018-2024 Zidan Zhang <zhangzidan@gmail.com>

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

import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
from collections import Counter


def getDomain(concentration, aP, aN, aPcoords, aNcoords, Dxyz):
    if concentration == 'c000':
        domains = [29.387, 9.635, 97.136, 10.542]  # c000
    elif concentration == 'c020':
        domains = [28.330, 11.155, 100.320, 10.950]  # c020
    elif concentration == 'c040':
        domains = [26.760, 10.590, 102.400, 11.142]  # c040
    elif concentration == 'c060':
        domains = [26.140, 10.776, 105.340, 10.446]  # c060
    elif concentration == 'c080':
        domains = [25.320, 10.667, 107.740, 11.429]  # c080
    if1 = [domains[0] - domains[1] / 2, domains[0] + domains[1] / 2]
    bk1 = [domains[0] + domains[1] / 2, domains[2] - domains[3] / 2]
    if2 = [domains[2] - domains[3] / 2, domains[2] + domains[3] / 2]
    to1 = [domains[0] - domains[1] / 2, domains[2] + domains[3] / 2]
    aPinif1 = []; aNinif1 = []
    aPinif2 = []; aNinif2 = []
    aPinbk  = []; aNinbk = []
    aPinto  = []; aNinto = []
    for i in range(len(aP)):
        if (aPcoords[:, Dxyz][i] >= if1[0] and aPcoords[:, Dxyz][i] <= if1[1]):
            aPinif1.append(i)
        if (aPcoords[:, Dxyz][i] >= if2[0] and aPcoords[:, Dxyz][i] <= if2[1]):
            aPinif2.append(i)
        if (aPcoords[:, Dxyz][i] >= bk1[0] and aPcoords[:, Dxyz][i] <= bk1[1]):
            aPinbk.append(i)
        if (aPcoords[:, Dxyz][i] >= to1[0] and aPcoords[:, Dxyz][i] <= to1[1]):
            aPinto.append(i)
    for i in range(len(aN)):
        if (aNcoords[:, Dxyz][i] >= if1[0] and aNcoords[:, Dxyz][i] <= if1[1]):
            aNinif1.append(i)
        if (aNcoords[:, Dxyz][i] >= if2[0] and aNcoords[:, Dxyz][i] <= if2[1]):
            aNinif2.append(i)
        if (aNcoords[:, Dxyz][i] >= bk1[0] and aNcoords[:, Dxyz][i] <= bk1[1]):
            aNinbk.append(i)
        if (aNcoords[:, Dxyz][i] >= to1[0] and aNcoords[:, Dxyz][i] <= to1[1]):
            aNinto.append(i)
    aPinif = aPinif1 + aPinif2
    aNinif = aNinif1 + aNinif2
    return aPinto, aPinbk, aPinif, aNinto, aNinbk, aNinif


def analysis(top, traj, anion, cation, cutoff, concentration, Dxyz, frame_id):
    uta = mda.Universe(top, traj)
    aP = uta.select_atoms("type " + anion)
    aN = uta.select_atoms("type " + cation) 
    ts = uta.trajectory[frame_id]
    cell = ts.dimensions
    distAC = dist.distance_array(aP.positions, aN.positions, cell)          
    htt_to, htt_bk, htt_if = [], [], []
    aPinto, aPinbk, aPinif, aNinto, aNinbk, aNinif = getDomain(concentration, aP, aN, aP.positions, aN.positions, Dxyz)
    for i in range(0, len(aP)):          
        for j in range(0, len(aN)):
            if distAC[i][j] < cutoff:
                if (i in aPinto) and (j in aNinto):
                    htt_to.append((i, j))
                if (i in aPinbk) and (j in aNinbk):
                    htt_bk.append((i, j))
                if (i in aPinif) and (j in aNinif):
                    htt_if.append((i, j))

    return ts.frame, htt_to, htt_bk, htt_if
