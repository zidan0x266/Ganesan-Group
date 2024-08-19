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
import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
import MDAnalysis.lib.distances as distances
import numpy as np
from analysis import *


"""
    assumptions:
    1. the microphase separated structure has a perfect/smooth interface
    2. particle distribution is homogeneous within specific domain
    3. rmax equals to the half of the specific domain
"""


def elementalV(resolution, rmax, bound, xz):
    """
    This function gets the elemental volume dV for the reference point
    resolution: number of bins
    rmax: cutoff for the rdf, could be half of the simulation box or half of the domain
    rinter: distance of the reference particle to the interface
    """
    nbins = np.int(resolution)
    rdfrange = (0.0, rmax)
    edges = np.histogram(rdfrange, nbins)[1]
    bins = 0.5 * (edges[:-1] + edges[1:])
    dr = edges[1]
    vol = np.power(edges[1:], 3) - np.power(edges[:-1], 3)
    vol *= 4/3.0 * np.pi
    if np.absolute(xz - bound[0]) < rmax:  # Spherical Cap
        rinter = np.absolute(xz - bound[0])
        upper = np.int(np.floor(rinter / dr))
        for i in range(upper, len(edges) - 1):
            vcap = 2 * np.pi * edges[i] * (edges[i] - rinter) * dr  # Scap = 2piRh
            vol[i] -= vcap
    if np.absolute(xz - bound[1]) < rmax:  # Spherical Cap
        linter = np.absolute(xz - bound[1])
        upper = np.int(np.floor(linter / dr))
        for i in range(upper, len(edges) - 1):
            vcap = 2 * np.pi * edges[i] * (edges[i] - linter) * dr  # Scap = 2piRh
            vol[i] -= vcap
    return vol


def rdfPrint(filename, bins, rdf):  # print displacement
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# r g(r)', file=anaout)
        for i in range(0, len(bins)):
            print('{:10.3f} {:10.5f}'.format(bins[i], rdf[i]), file=anaout)

            
uta = mda.Universe('../cg_topol.tpr', '../cg_traj.xtc')
aP = uta.select_atoms("type PF")
aN = uta.select_atoms("type IM")
#if1 = [27.5, 32.5]
#bk1 = [32.5, 47.5]
#if2 = [47.5, 52.5]
#to1 = [27.5, 52.5]
if1 = [52.5, 57.5]
bk1 = [57.5, 72.5]
if2 = [72.5, 77.5]
to1 = [52.5, 77.5]
rdf1 = []
rdf2 = []
rdf3 = []
rdf4 = []
for ts in uta.trajectory:
    if ts.frame >= 10000:
        cell = uta.trajectory.ts.dimensions
        Ainif1 = []; Cinif1 = []
        Ainif2 = []; Cinif2 = []
        Ainbk = []; Cinbk = []
        Ainto = []; Cinto = []
        for i in range(len(aP)):
            if (aP.positions[:,2][i] >= if1[0] and aP.positions[:,2][i] <= if1[1]):
                Ainif1.append(i)
            if (aP.positions[:,2][i] >= if2[0] and aP.positions[:,2][i] <= if2[1]):
                Ainif2.append(i)
            if (aP.positions[:,2][i] >= bk1[0] and aP.positions[:,2][i] <= bk1[1]):
                Ainbk.append(i)
            if (aP.positions[:,2][i] >= to1[0] and aP.positions[:,2][i] <= to1[1]):
                Ainto.append(i)
            if (aN.positions[:,2][i] >= if1[0] and aN.positions[:,2][i] <= if1[1]):
                Cinif1.append(i)
            if (aN.positions[:,2][i] >= if2[0] and aN.positions[:,2][i] <= if2[1]):
                Cinif2.append(i)
            if (aN.positions[:,2][i] >= bk1[0] and aN.positions[:,2][i] <= bk1[1]):
                Cinbk.append(i)
            if (aN.positions[:,2][i] >= to1[0] and aN.positions[:,2][i] <= to1[1]):
                Cinto.append(i)
        volume_if = cell[0] ** 2 * 5
        volume_bk = cell[0] ** 2 * 15
        volume_to = cell[0] ** 2 * 25
        den_if1 = len(Ainif1) * len(Cinif1) / volume_if
        den_if2 = len(Ainif2) * len(Cinif2) / volume_if
        den_bk = len(Ainbk) * len(Cinbk) / volume_bk
        den_to = len(Ainto) * len(Cinto) / volume_to
        aNposit1 = np.empty((len(Cinif1), 3))
        aNposit2 = np.empty((len(Cinif2), 3))
        aNposit3 = np.empty((len(Cinbk), 3))
        aNposit4 = np.empty((len(Cinto), 3))
        i = 0; j = 0; k = 0; l = 0
        for aid in Cinif1:
            aNposit1[i] = aN.positions[aid]
            i += 1
        for aid in Cinif2:
            aNposit2[j] = aN.positions[aid]
            j += 1
        for aid in Cinbk:
            aNposit3[k] = aN.positions[aid]
            k += 1
        for aid in Cinto:
            aNposit4[l] = aN.positions[aid]
            l += 1
        tmpif1 = []
        for i in Ainif1:
            distpair1 = distances.capped_distance(aP.positions[i], aNposit1, 15, box = cell)[1]
            count = np.histogram(distpair1, 200, (0, 15))[0]
            vol = elementalV(200, 15, if1, aP.positions[i][2])
            tmpif1.append(count / (den_if1 * vol))
        tmpif2 = []
        for i in Ainif2:
            distpair2 = distances.capped_distance(aP.positions[i], aNposit2, 15, box = cell)[1]
            count = np.histogram(distpair2, 200, (0, 15))[0]
            vol = elementalV(200, 15, if2, aP.positions[i][2])
            tmpif2.append(count / (den_if2 * vol))
        tmpbk = []
        for i in Ainbk:
            distpair3 = distances.capped_distance(aP.positions[i], aNposit3, 15, box = cell)[1]
            count = np.histogram(distpair3, 200, (0, 15))[0]
            vol = elementalV(200, 15, bk1, aP.positions[i][2])
            tmpbk.append(count / (den_bk * vol))
        tmpto = []
        for i in Ainto:
            distpair4 = distances.capped_distance(aP.positions[i], aNposit4, 15, box = cell)[1]
            count = np.histogram(distpair4, 200, (0, 15))[0]
            vol = elementalV(200, 15, to1, aP.positions[i][2])
            tmpto.append(count / (den_to * vol))
        rdf1.append(np.sum(tmpif1, axis = 0))
        rdf2.append(np.sum(tmpif2, axis = 0))
        rdf3.append(np.sum(tmpbk, axis = 0))
        rdf4.append(np.sum(tmpto, axis = 0))
    if ts.frame == 15000:
        break
rdf1 = np.average(rdf1, axis = 0)
rdf2 = np.average(rdf2, axis = 0)
rdf_if = np.average([rdf1, rdf2], axis = 0)
rdf_bk = np.average(rdf3, axis = 0)
rdf_to = np.average(rdf4, axis = 0)
nbins = np.int(200)
rdfrange = (0.0, 15)
edges = np.histogram(rdfrange, nbins)[1]
bins = 0.5 * (edges[:-1] + edges[1:])
rdfPrint('nonc_if', bins, rdf_if)
rdfPrint('nonc_bk', bins, rdf_bk)
rdfPrint('nonc_to', bins, rdf_to)