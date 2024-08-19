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

import multiprocessing as mp
import functools

import MDAnalysis as mda
import numpy as np


def removePBC(cell, curCoords, prevCoords):
    unwrapCoords = curCoords.copy()
    dim = [cell[i] for i in range(0, 3)]
    hbox = [0.5 * dim[i] for i in range(0, 3)]
    natoms = len(curCoords)
    shiftx = np.zeros(natoms)
    shifty = np.zeros(natoms)
    shiftz = np.zeros(natoms)
    tmpdz = curCoords[:,2] - prevCoords[:,2]  # only consider z-direction at this time
    i = 0
    for dz in tmpdz:
        shift = 0
        if dz > hbox[2]:
            shift -= 1
            while (dz + shift * dim[2]) > hbox[2]:
                shift -= 1
        elif dz < -hbox[2]:
            shift += 1
            while (dz + shift * dim[2]) < -hbox[2]:
                shift += 1
        shiftz[i] = shift
        i += 1
    unwrapCoords[:,2] += (shiftz * dim[2])
    tmpdy = curCoords[:,1] - prevCoords[:,1]  # only consider y-direction at this time
    i = 0
    for dy in tmpdy:
        shift = 0
        if dy > hbox[1]:
            shift -= 1
            while (dy + shift * dim[1]) > hbox[1]:
                shift -= 1
        elif dy < -hbox[1]:
            shift += 1
            while (dy + shift * dim[1]) < -hbox[1]:
                shift += 1
        shifty[i] = shift
        i += 1
    unwrapCoords[:,1] += (shifty * dim[1])
    tmpdx = curCoords[:,0] - prevCoords[:,0]  # only consider x-direction at this time
    i = 0
    for dx in tmpdx:
        shift = 0
        if dx > hbox[0]:
            shift -= 1
            while (dx + shift * dim[0]) > hbox[0]:
                shift -= 1
        elif dx < -hbox[0]:
            shift += 1
            while (dx + shift * dim[0]) < -hbox[0]:
                shift += 1
        shiftx[i] = shift
        i += 1
    unwrapCoords[:,0] += (shiftx * dim[0])        
    return unwrapCoords


def main():
    uta = mda.Universe("md.tpr", "test0.xtc")
    frames = uta.trajectory
    unwrap = mda.Writer("ps100.xtc", frames.n_atoms)
    firstframe = True
    for ts in uta.trajectory:
        cell = ts.dimensions
        if firstframe:
            prevaL = ts.positions
            firstframe = False
            unwrap.write(uta.atoms)
        else:
            curaL = ts.positions
            tmpaL = removePBC(cell, curaL, prevaL)
            prevaL = tmpaL
            ts._pos[:] = tmpaL[:]
            if ts.frame % 100 == 0:
                unwrap.write(uta.atoms)


if __name__ == "__main__":
    main()
