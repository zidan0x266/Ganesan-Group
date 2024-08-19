"""
Copyright (C) 2018-2019 Zidan Zhang <zhangzidan@gmail.com>
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
import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
from collections import Counter


def analysis(top, traj, anion, cation, beadiny, cutoff, Icutoff, counterN, frame_id):
    """
    This function outputs the analysis for everything of interest for
    copolymer (random or block), including association relationship,
    time auto correlation function and these properties near the
    interface.

    Process a single frame.
    """
    uta = mda.Universe(top, traj)  
    # Group atoms to be analyzied
    aP = uta.select_atoms("type " + anion)
    aN = uta.select_atoms("type " + cation)
    aE = uta.select_atoms("type " + beadiny)  
    ts = uta.trajectory[frame_id]
    cell = uta.trajectory.ts.dimensions 
    assoAC = []; assoCH = []; assoAE = [] # total associations
    assoBAC = []; assoBCH = []
    assoIAC = []; assoIAE = []; assoICH = []
    distAC = dist.distance_array(aP.positions, aN.positions, cell)   
    distAE = dist.distance_array(aP.positions, aE.positions, cell)           
    htt = [] # Pair events at time now
    sub_htt = [] # Pair events at time now near interface
    for i in range(0, len(aP)):
        beads = []; chains = []
        t1 = 0; t2 = 0               
        for j in range(0, len(aN)):
            if distAC[i][j] < cutoff:
                t1 += 1
                beads.append(j)
                t4 = (i, j)
                htt.append(t4)
                for k in range(0, len(aE)):
                    if distAE[i][k] < Icutoff:
                        t5 = (i, j)
                        sub_htt.append(t5)
                        t2 += 1                 
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
        assoAE.append(t2)
        if assoAC[i] != 0 and assoAE[i] != 0:
            assoIAE.append(i)
    for i in assoIAE:
        assoIAC.append(assoAC[i])
        assoICH.append(assoCH[i])
    taP = np.linspace(0, len(aP) - 1, len(aP), dtype = int)
    aPinBulk = [i for i in taP if i not in assoIAE]
    for i in aPinBulk:
        assoBAC.append(assoAC[i])
        assoBCH.append(assoCH[i])
    it_htt = set(sub_htt)
    in_htt = set(htt) - it_htt
    return ts.frame, Counter(assoAC), Counter(assoCH), Counter(assoBAC), Counter(assoBCH), Counter(assoIAC), Counter(assoICH), len(assoIAE), htt, list(in_htt), list(it_htt)
