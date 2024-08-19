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
    top = "../cg_topol.tpr"
    trj = "../mergextc/subdata.xtc"
    uta = mda.Universe(top, trj)
    aT1 = uta.select_atoms("type PF")
    aT2 = uta.select_atoms("type IM")
    rdf1 = ionpair.InterRDF(aT1, aT1, nbins=1000, range=(0.0, 20.0), exclusion_block=(1,1))
    rdf1.run()
    rdfPrint('cgrdf_PF_PF', rdf1.bins, rdf1.rdf)
    rdf2 = ionpair.InterRDF(aT2, aT2, nbins=1000, range=(0.0, 20.0), exclusion_block=(1,1))
    rdf2.run()
    rdfPrint('cgrdf_IM_IM', rdf2.bins, rdf2.rdf)
    rdf3 = ionpair.InterRDF(aT1, aT2, nbins=1000, range=(0.0, 20.0), exclusion_block=(1,1))
    rdf3.run()
    rdfPrint('cgrdf_PF_IM', rdf3.bins, rdf3.rdf)


if __name__ == "__main__":
    main()
