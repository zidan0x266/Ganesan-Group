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

def dpPrint(filename, dt, dx):  # print displacement
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# time dx', file=anaout)
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

uta = mda.Universe("cg_topol.tpr", "cg_unwrap.xtc")

aL = uta.select_atoms('all')  # all atoms

al_cod = []
for ts in uta.trajectory:
    if ts.frame >= 5000:
        al_cod.append(aL.center_of_mass()[0])

dt = np.linspace(0, 15000, 15001)
dx_al = finalGetdx_i(al_cod, 0, 1, 20)
dpPrint('dx_aL', dt, dx_al)

print("Displacement analysis finished!")