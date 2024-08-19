"""
Copyright (C) 2018-2021 Zidan Zhang <zhangzidan@gmail.com>
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
import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
from collections import Counter


def analysis(top, traj, anion, cation, cutoff, counterN, frame_id):
    """
    This function outputs the analysis for everything of interest for
    homo polymer, including association relationship,
    time auto correlation function.

    Process a single frame.
    """
    uta = mda.Universe(top, traj, format="LAMMPS")
    aP = uta.select_atoms("type " + anion)
    aN = uta.select_atoms("type " + cation) 
    ts = uta.trajectory[frame_id]
    cell = ts.dimensions
    assoAC = []; assoCH = [] # total associations
    distAC = dist.distance_array(aP.positions, aN.positions, cell)          
    htt = [] # Pair events at time now
    for i in range(0, len(aP)):
        beads = []; chains = []
        t1 = 0              
        for j in range(0, len(aN)):
            if distAC[i][j] < cutoff:
                t1 += 1
                beads.append(j)
                htt.append((i, j))
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
    return ts.frame, Counter(assoAC), Counter(assoCH), htt
