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


def dpPrint(filename, dt, dx):  # print displacement
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# time delta-X', file=anaout)
        for i in range(0, len(dt)):
            print('{:10.3f} {:10.5f}'.format(dt[i], dx[i]), file=anaout)


def getdx_i(COORDs, t0, t):
    """
    This function give the intermitent Association Lifetime
    C_i = <h(t0)h(t)>/<h(t0)> between t0 and t
    """
    C_i = COORDs[t] - COORDs[t0]
    return C_i


def intervdx_i(COORDs, t0, tf, reservoir, dt):
    """
    This function gets all the data for the h(t0)h(t0+dt), where
    t0 = 1,2,3,...,tf. This function give us a point of the final plot
    C(t) vs t,
    """
    a = 0
    count = 0
    for t0 in reservoir:
        if t0 + dt <= tf:
            a += getdx_i(COORDs, t0, t0 + dt)
            count += 1
        else:
            break
    return a / count


def finalGetdx_i(COORDs, t0, trestart, nt):
    """
    This function gets the final data of the C_i graph.
    """
    tf = len(COORDs) - 1
    num_points = len(COORDs)
    pool = mp.Pool(processes=nt)
    reservoir = [x * trestart for x in range(tf) if x * trestart < tf]
    func = functools.partial(intervdx_i, COORDs, t0, tf, reservoir)
    return_data = pool.map(func, list(range(num_points)))
    return return_data


uta = mda.Universe("cg_topol.tpr", "cg_traj.xtc")

aL = uta.select_atoms('all')  # all atoms
aP = uta.select_atoms('type PF')  # anion PF6
aN = uta.select_atoms('type IM')  # anion PF6

axis = 0
al_unwrap = []
ap_unwrap = []
an_unwrap = []
firstframe = True
for ts in uta.trajectory:
    if ts.frame >= 10000:
        cell = ts.dimensions
        if firstframe:
            prevaL = aL.positions
            prevaP = aP.positions
            prevaN = aN.positions
            firstframe = False
        curaL = aL.positions
        curaP = aP.positions
        curaN = aN.positions
        tmpaL = removePBC(cell, curaL, prevaL)
        tmpaP = removePBC(cell, curaP, prevaP)
        tmpaN = removePBC(cell, curaN, prevaN)
        aL.positions = tmpaL
        com = aL.centroid()
        prevaP = tmpaP
        prevaN = tmpaN
        ap_unwrap.append(tmpaP[:,axis] - com[axis])
        an_unwrap.append(tmpaN[:,axis] - com[axis])
        #print("Frame {}".format(ts.frame))
    if ts.frame % 1000 == 0:
        print("Frame {}".format(ts.frame))
        #break

dt = np.linspace(0, 30000, 30001)
dxap = finalGetdx_i(ap_unwrap, 0, 1, 24)
dx_ap = np.average(dxap, axis = 1)
dpPrint('dx_nocom_aP', dt, dx_ap)
dxan = finalGetdx_i(an_unwrap, 0, 1, 24)
dx_an = np.average(dxan, axis = 1)
dpPrint('dx_nocom_aN', dt, dx_an)

print("Displacement analysis finished!")

