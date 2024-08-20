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

import multiprocessing as mp
import functools
import numpy as np
import h5py
import time
import MDAnalysis as mda


def travelPrint(filename, distance):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# distance', file=anaout)
        for i in range(len(distance)):
            print('{:10.5f}'.format(distance[i]), file=anaout)


def main():
    # generic parameters
    nframes = 15001  # number of frames
    nt = 56  # number of processors
    tstar = 1147  # distance diff shift
    # generate the coordinates object
    top = "../cg_topol.tpr"
    trj = "../fong/unwrap.xtc"
    atp = "PF"  # fetch central particle
    uta = mda.Universe(top, trj)
    atg = uta.select_atoms("type " + atp)
    COORDs = []
    for ts in uta.trajectory:
        if ts.frame % 1000 == 0:
            print("The current reference frame is: " + str(ts.frame))
        COORDs.append(atg.positions)
        if ts.frame == nframes:
            break
    START = np.linspace(0, 10000, 6, dtype = int)
    travel_length = []
    for start in START:
        print("The current reference frame is: " + str(start))
        if (start + tstar) <= nframes:
            displ = COORDs[start + tstar] - COORDs[start]
            displacement = np.sqrt(np.sum(np.square(displ), axis = 1))
        for length in displacement:
            travel_length.append(length)
    travelPrint('length_till_tstar', travel_length)

if __name__ == "__main__":
    main()
