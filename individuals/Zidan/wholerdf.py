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


"""
    assumptions:
    1. the microphase separated structure has a perfect/smooth interface
    2. particle distribution is homogeneous within specific domain
    3. rmax equals to the half of the specific domain
"""


def elementalV(resolution, rmax):
    """
    This function gets the elemental volume dV for the reference point
    resolution: number of bins
    rmax: cutoff for the rdf, could be half of the simulation box or half of the domain
    """
    nbins = np.int(resolution)
    rdfrange = (0.0, rmax)
    edges = np.histogram(rdfrange, nbins)[1]
    bins = 0.5 * (edges[:-1] + edges[1:])
    dr = edges[1]
    vol = np.power(edges[1:], 3) - np.power(edges[:-1], 3)
    vol *= 4/3.0 * np.pi

    return vol


def rdfPrint(filename, bins, rdf):  # print displacement
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# r g(r)', file=anaout)
        for i in range(0, len(bins)):
            print('{:10.3f} {:10.5f}'.format(bins[i], rdf[i]), file=anaout)


uta = mda.Universe('../cg_topol.tpr', '../cg_traj.xtc')
aP = uta.select_atoms("type PF")
aN = uta.select_atoms("type IM")
N = len(aP) * len(aN)
volume = uta.trajectory.ts.volume # volume fraction
density = N / volume
rdfs = []
for ts in uta.trajectory:
    if ts.frame >= 60000:
        cell = uta.trajectory.ts.dimensions
        tmp = []
        for i in range(len(aP)):
            distpair = distances.capped_distance(aP.positions[i], aN.positions, 15, box = cell)[1]
            count = np.histogram(distpair, 200, (0, 15))[0]
            #vol = elementalV(resolution, 15)
            #tmp.append(count / (density * vol))
            tmp.append(count)
        rdfs.append(np.average(tmp, axis = 0))
    if ts.frame == 60100:
        break

rdf = np.average(rdfs, axis = 0)
cum_rdf = np.cumsum(rdf)
nbins = np.int(200)
rdfrange = (0.0, 15)
edges = np.histogram(rdfrange, nbins)[1]
bins = 0.5 * (edges[:-1] + edges[1:])
#rdfPrint('rdf', bins, rdf)
rdfPrint('rdf_ho', bins, cum_rdf)