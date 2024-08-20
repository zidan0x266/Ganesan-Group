"""
Copyright (C) 2018-2020 Zidan Zhang <zhangzidan@gmail.com>
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
import numpy as np


def dihPrint(filename, acf):  # print sqt to the output
    time = np.linspace(0, len(acf) - 1, len(acf))
    with open('{}.xvg'.format(filename), 'w') as sqtout:
        print('# time, dihedral auto-correlation function', file=sqtout)
        for i in range(len(time)):
            print('{:10.3f} {:10.5f}'.format(time[i], acf[i]), file=sqtout)


def main():
    uta = mda.Universe('md.tpr', 'ps100.xtc')
    backbone = uta.select_atoms('name C1 or name C2')  # select atoms only on backbone
    dihindex = []
    for i, dih in enumerate(backbone.dihedrals):
        if dih.atoms[0].name == 'C1' and dih.atoms[1].name == 'C2' and dih.atoms[2].name == 'C1' and dih.atoms[-1].name == 'C2':
            dihindex.append(i)
        #if dih.atoms[0].name == 'C2' and dih.atoms[1].name == 'C1' and dih.atoms[2].name == 'C2' and dih.atoms[-1].name == 'C1':
        #    dihindex.append(i)
    nframes = len(uta.trajectory)
    cosvals = np.zeros((nframes, len(dihindex)))
    acfphi = np.zeros(nframes)
    firstframe = True
    for ts in uta.trajectory:
        frame = ts.frame
        if frame % 10 == 0:
            print('Now processing percentage: {:5.3f}'.format(float(frame) / (nframes - 1)))
        for i in range(len(dihindex)):
            phi = backbone.dihedrals.values()[dihindex[i]]
            cosvals[frame][i] = np.cos(phi)
        if firstframe:
            avg1 = np.average([i ** 2 for i in cosvals[frame]])
            avg2 = np.average(cosvals[frame]) ** 2
            firstframe = False
        avg3 = np.average([i * j for i, j in zip(cosvals[frame], cosvals[0])])
        acfphi[frame] = (avg3 - avg2) / (avg1 - avg2)
    dihPrint('acf_dih', acfphi)
        

if __name__ == "__main__":
    main()