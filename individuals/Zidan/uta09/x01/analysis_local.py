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


def analysis(top, traj, anion, cation, cutoff, frame_id):
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
    # Ranges for local definition
    #ifm1 = [26.9, 31.9]; ifn1 = [71.9, 76.9]
    #bkm1 = [31.9, 48.1]; bkn1 = [76.9, 93.1]
    #ifm2 = [48.1, 53.1]; ifn2 = [93.1, 98.1]
    #tom1 = [26.9, 53.1]; ton1 = [71.9, 98.1]
    ifm1 = [8.0]; ifn1 = [8.0, 26.9]
    bkm1 = [26.9, 53.1]; bkn1 = [76.9, 100.1]
    ifm2 = [53.1, 76.9]; ifn2 = [100.1]
    # Load frame
    ts = uta.trajectory[frame_id]
    cell = uta.trajectory.ts.dimensions
    # local range processing
    Ainif1 = []; Ainif2 = []; Ainif3 = []; Ainif4 = []  # Atoms in interface
    Ainbk1 = []; Ainbk2 = []  # Atoms in bulk
    for i in range(len(aP)):
        # interface
        if (aP.positions[:,2][i] < ifm1[0]):
            Ainif1.append(i)
        if (aP.positions[:,2][i] >= ifm2[0] and aP.positions[:,2][i] < ifm2[1]):
            Ainif2.append(i)
        if (aP.positions[:,2][i] >= ifn1[0] and aP.positions[:,2][i] < ifn1[1]):
            Ainif3.append(i)
        if (aP.positions[:,2][i] >= ifn2[0]):
            Ainif4.append(i)
        # bulk
        if (aP.positions[:,2][i] >= bkm1[0] and aP.positions[:,2][i] < bkm1[1]):
            Ainbk1.append(i)
        if (aP.positions[:,2][i] >= bkn1[0] and aP.positions[:,2][i] < bkn1[1]):
            Ainbk2.append(i)
    Ainif = sorted(Ainif1 + Ainif2 + Ainif3 + Ainif4)
    Ainbk = sorted(Ainbk1 + Ainbk2)
    Ainto = sorted(Ainif + Ainbk)
    Zinif1 = []; Zinif2 = []; Zinif3 = []; Zinif4 = []  # Atoms in interface
    Zinbk1 = []; Zinbk2 = []  # Atoms in bulk
    for i in range(len(aN)):
        # interface
        if (aN.positions[:,2][i] < ifm1[0]):
            Zinif1.append(i)
        if (aN.positions[:,2][i] >= ifm2[0] and aN.positions[:,2][i] < ifm2[1]):
            Zinif2.append(i)
        if (aN.positions[:,2][i] >= ifn1[0] and aN.positions[:,2][i] < ifn1[1]):
            Zinif3.append(i)
        if (aN.positions[:,2][i] >= ifn2[0]):
            Zinif4.append(i)
        # bulk
        if (aN.positions[:,2][i] >= bkm1[0] and aN.positions[:,2][i] < bkm1[1]):
            Zinbk1.append(i)
        if (aN.positions[:,2][i] >= bkn1[0] and aN.positions[:,2][i] < bkn1[1]):
            Zinbk2.append(i)
    Zinif = sorted(Zinif1 + Zinif2 + Zinif3 + Zinif4)
    Zinbk = sorted(Zinbk1 + Zinbk2)
    Zinto = sorted(Zinif + Zinbk)
    # Postprocessing
    distAC = dist.distance_array(aP.positions, aN.positions, cell)           
    htt = [] # Pair events at time now
    if_htt = [] # Pair events at time now near interface
    bk_htt = [] # Pair events at time now in bulk
    to_htt = [] # Pair events at time now in AA-riched phase
    for i in range(0, len(aP)):
        beads = []
        t1 = 0              
        for j in range(0, len(aN)):
            if distAC[i][j] <= cutoff:
                t1 += 1
                beads.append(j)
                t4 = (i, j)
                htt.append(t4)
                if i in Ainif and j in Zinif:
                    if_htt.append(t4)
                if i in Ainbk and j in Zinbk:
                    bk_htt.append(t4)
                if i in Ainto and j in Zinto:
                    to_htt.append(t4)                  
    return ts.frame, htt, if_htt, bk_htt, to_htt

