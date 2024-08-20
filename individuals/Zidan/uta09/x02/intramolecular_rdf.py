"""
Copyright (C) 2018-2022 Zidan Zhang <zhangzidan@gmail.com>

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


def main():
    top = '../cg_topol.tpr'
    trj = '../cgtrj.xtc'
    ref = "TF"
    sel = "IM"
    uta = mda.Universe(top, trj)
    aP = uta.select_atoms("type " + ref)
    aN = uta.select_atoms("type " + sel)
    rdfmax = 20
    resolution = 200
    N = len(aP) * len(aN)
    volume = uta.trajectory.ts.volume # volume fraction
    density = N / volume
    rdfs = []
    for ts in uta.trajectory:
        cell = uta.trajectory.ts.dimensions
        tmp = []
        for i in range(len(aP)):
            distpair = distances.capped_distance(aP.positions[i], aN.positions, rdfmax, box = cell)[1]
            count = np.histogram(distpair, resolution, (0, rdfmax))[0]
            vol = elementalV(resolution, rdfmax)
            tmp.append(count / (density * vol))
        rdfs.append(np.average(tmp, axis = 0))
    rdf = np.average(rdfs, axis = 0)
    nbins = np.int(resolution)
    rdfrange = (0.0, rdfmax)
    edges = np.histogram(rdfrange, nbins)[1]
    bins = 0.5 * (edges[:-1] + edges[1:])
    rdfPrint('rdf', bins, rdf)


if __name__ == "__main__":
    main()