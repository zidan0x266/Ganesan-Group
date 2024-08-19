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

import multiprocessing as mp
import functools
import numpy as np
import MDAnalysis as mda

def main():
    top = '../cg_topol.tpr'
    trj = '../cgtrj.xtc'
    ref = "TF"
    Dxyz = 2  # 0: X; 1: Y; 2: Z
    cases = ['LA01', 'LA02', 'LA03']
    curr_case = cases[2]
    if curr_case == cases[0]:
        domains = [14.18, 32.86, 66.80, 87.46]  # A30B10
    elif curr_case == cases[1]:
        domains = [22.06, 32.38, 83.64, 92.90]  # A45B15
    elif curr_case == cases[2]:
        domains = [25.42, 33.96, 93.06, 100.94]  # A60B20
    uta = mda.Universe(top, trj)
    atg = uta.select_atoms("type " + ref)
    a1, a2= [], []
    for ts in uta.trajectory:
        t1 = len([atom for atom in range(len(atg)) if atg.positions[atom][Dxyz] >= domains[0] and atg.positions[atom][Dxyz] < domains[1]])
        t2 = len([atom for atom in range(len(atg)) if atg.positions[atom][Dxyz] >= domains[1] and atg.positions[atom][Dxyz] < domains[2]])
        t3 = len([atom for atom in range(len(atg)) if atg.positions[atom][Dxyz] >= domains[2] and atg.positions[atom][Dxyz] < domains[3]])
        a1.append(t2 / len(atg))  # fraction of ions in the bulk region
        a2.append((t1 + t3) / len(atg))  # fraction of ions in the interfacial region
    print("{:8.5f}".format(np.average(a1)))
    print("{:8.5f}".format(np.std(a1)))
    print("{:8.5f}".format(np.average(a2)))
    print("{:8.5f}".format(np.std(a2)))


if __name__ == "__main__":
    main()