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
import h5py
import time
import MDAnalysis as mda

def main():
    range00 = [2.672042, 3.719957, 5.625372, 7.530788, 9.436204, 10.475795]
    range02 = [2.649564, 3.618435, 5.669264, 7.720093, 9.770923, 10.871076]
    range04 = [2.516654, 3.529345, 5.683280, 7.837216, 9.991152, 11.080847]
    range06 = [2.361888, 3.548111, 5.793017, 8.037923, 10.282829, 11.369170]
    range08 = [2.422559, 3.415440, 5.751478, 8.087515, 10.423554, 11.616445]
    domains = [x * 10 for x in range08]
    top = "cg_topol.tpr"
    trj = "data.xtc"
    uta = mda.Universe(top, trj)
    atp = "TF"
    atg = uta.select_atoms("type " + atp)
    axis = 1  # 0: X, 1: Y; 2: Z
    a1, a2, a3, a4, a5 = [], [], [], [], []
    for ts in uta.trajectory:
        t1, t2, t3, t4, t5 = 0, 0, 0, 0, 0
        for atom in range(len(atg)):
            if domains[0] <= atg.positions[atom][axis] < domains[1]:
                t1 += 1
            elif domains[1] <= atg.positions[atom][axis] < domains[2]:
                t2 += 1
            if domains[2] <= atg.positions[atom][axis] < domains[3]:
                t3 += 1
            if domains[3] <= atg.positions[atom][axis] < domains[4]:
                t4 += 1
            if domains[4] <= atg.positions[atom][axis] < domains[5]:
                t5 += 1
        a1.append(t1)
        a2.append(t2)
        a3.append(t3)
        a4.append(t4)
        a5.append(t5)
    print("{:8.3f}, {:8.3f}, {:8.3f}, {:8.3f}, {:8.3f}".format(np.average(a1), np.average(a2), np.average(a3), np.average(a4), np.average(a5)))
    print("{:8.3f}, {:8.3f}, {:8.3f}, {:8.3f}, {:8.3f}".format(np.std(a1), np.std(a2), np.std(a3), np.std(a4), np.std(a5)))


if __name__ == "__main__":
    main()