"""
Copyright (C) 2018-2021 Zidan Zhang <zhangzidan@gmail.com>

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
    top = "../md.tpr"
    trj = "../md.trr"
    start = 30000
    uta = mda.Universe(top, trj)
    aR = uta.select_atoms("type cu")  # atoms in the reference group
    aS1 = uta.select_atoms("resname MPIG and name O3")  # amide oxygen
    aS2 = uta.select_atoms("resname TFSI and name NA1")  # TFSI nitrogen
    aS3 = uta.select_atoms("resname MPIG and name N3")  # imidazole nitrogen
    aS4 = uta.select_atoms("resname HPEO and name O1 or resname MPEO and name O1 or resname TPEO and name O1")  # PEO oxygen
    rdfs = []
    for ts in uta.trajectory:
        if ts.frame >= start:
            cell = uta.trajectory.ts.dimensions
            tmp = []
            for i in range(len(aR)):
                distpair = distances.capped_distance(aR.positions[i], aS4.positions, 15, box = cell)[1]
                count = np.histogram(distpair, 200, (0, 15))[0]
                tmp.append(count)
            rdfs.append(np.average(tmp, axis = 0))
        if ts.frame == start + 100:
            break
    
    rdf = np.average(rdfs, axis = 0)
    cum_rdf = np.cumsum(rdf)
    nbins = int(200)
    rdfrange = (0.0, 15)
    edges = np.histogram(rdfrange, nbins)[1]
    bins = 0.5 * (edges[:-1] + edges[1:])
    rdfPrint('cn_aS4', bins, cum_rdf)


if __name__ == "__main__":
    main()