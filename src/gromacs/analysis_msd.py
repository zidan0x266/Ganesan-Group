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


def msdPrint(filename, time, msd):  # print displacement
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# time delta-X', file=anaout)
        for i in range(0, len(time)):
            print('{:10.3f} {:10.5f}'.format(time[i], msd[i]), file=anaout)


def getdx_i(COORDs, t0, t):
    displ = COORDs[t] - COORDs[t0]
    msd = np.sum(np.square(displ), axis = 1)
    return np.average(msd)


def intervdx_i(COORDs, tf, reservoir, dt):
    a = 0
    count = 0
    for t0 in reservoir:
        if t0 + dt <= tf:
            a += getdx_i(COORDs, t0, t0 + dt)
            count += 1
        else:
            break
    return a / count


def finalGetdx_i(COORDs, trestart, nt):
    """
    This function gets the final data of the C_i graph.
    """
    tf = len(COORDs) - 1
    num_points = len(COORDs)
    pool = mp.Pool(processes=nt)
    reservoir = [x * trestart for x in range(tf) if x * trestart < tf]
    func = functools.partial(intervdx_i, COORDs, tf, reservoir)
    return_data = pool.map(func, list(range(num_points)))
    return return_data


uta = mda.Universe("../md.tpr", "unwrap.xtc")

aP = uta.select_atoms('type opls_480')  # anion PF6
COORDs = []
for ts in uta.trajectory:
    if ts.frame >= 0:  # start from arbitrary given reference frame
        COORDs.append(aP.positions)

dt = np.linspace(0, 10000, 10001)
dxap = finalGetdx_i(COORDs, 10, 2)
msdPrint('test', dt, dxap)

print("MSD analysis finished!")

