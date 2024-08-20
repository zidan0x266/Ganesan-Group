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
    if1 = [20.0, 30.0]; if2 = [47.0, 60.0]; if3 = [68.0, 80.0]; if4 = [97.0]
    bk1 = [30.0, 47.0]; bk2 = [80.0, 97.0]
    # Load frame
    ts = uta.trajectory[frame_id]
    cell = uta.trajectory.ts.dimensions
    # local range processing
    Ainif1 = []; Ainif2 = []; Ainif3 = []; Ainif4 = []  # Atoms in interface
    Ainbk1 = []; Ainbk2 = []  # Atoms in bulk
    for i in range(len(aP)):
        # interface
        if (aP.positions[:,2][i] >= if1[0] and aP.positions[:,2][i] < if1[1]):
            Ainif1.append(i)
        if (aP.positions[:,2][i] >= if2[0] and aP.positions[:,2][i] < if2[1]):
            Ainif2.append(i)
        if (aP.positions[:,2][i] >= if3[0] and aP.positions[:,2][i] < if3[1]):
            Ainif3.append(i)
        if (aP.positions[:,2][i] >= if4[0]):
            Ainif4.append(i)
        # bulk
        if (aP.positions[:,2][i] >= bk1[0] and aP.positions[:,2][i] < bk1[1]):
            Ainbk1.append(i)
        if (aP.positions[:,2][i] >= bk2[0] and aP.positions[:,2][i] < bk2[1]):
            Ainbk2.append(i)
    Ainif = sorted(Ainif1 + Ainif2 + Ainif3 + Ainif4)
    Ainbk = sorted(Ainbk1 + Ainbk2)
    Ainto = sorted(Ainif + Ainbk)
    # Postprocessing
    distAC = dist.distance_array(aP.positions, aN.positions, cell)           
    assoAC = []; assoSto = []; assoSbk = []; assoSif = []
    for i in range(0, len(aP)):
        t1 = 0              
        for j in range(0, len(aN)):
            if distAC[i][j] < cutoff:
                t1 += 1
        assoAC.append(t1)
        if i in Ainto:
            assoSto.append(t1)
        if i in Ainbk:
            assoSbk.append(t1)
        if i in Ainif:
            assoSif.append(t1)                 
    return ts.frame, Counter(assoAC), Counter(assoSto), Counter(assoSbk), Counter(assoSif)

