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

import MDAnalysis as mda
import MDAnalysis.analysis.rdf as ionpair

def rdfPrint(filename, bins, rdf):  # print RDF
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# r g(r)', file=anaout)
        for i in range(0, len(bins)):
            print('{:10.3f} {:10.5f}'.format(bins[i] / 10.0, rdf[i]), file=anaout)


def main():
    top = "../initial.data"
    trj = "../peo_stat_pige02.dcd"
    uta = mda.Universe(top, trj, format="LAMMPS")
    aR = uta.select_atoms("type 59")  # atoms in the reference group
    aS = uta.select_atoms("type 39")  # atoms in the selection group
    # type 59: lithium
    # type 39: amide oxygen
    # type 51: imidazolium nitrogen
    rdf = ionpair.InterRDF(aR, aS, nbins=1000, range=(0.0, 20.0))
    rdf.run()
    rdfPrint('rdf_lin3', rdf.bins, rdf.rdf)


if __name__ == "__main__":
    main()