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


def intervdx_i(COORDs, t0, tf, reservoir, dt):
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


def main():
    # Open Simulation
    top = "../initial.data"
    trj = "unwrap.dcd"
    uta = mda.Universe(top, trj, format="LAMMPS")
    aL = uta.select_atoms('all')  # all atom
    Li = uta.select_atoms('type 59')  # lithium
    TFSI = uta.select_atoms('type 54')  # lithium
    COORDs_Li = []
    COORDs_TFSI = []
    for ts in uta.trajectory:
        if ts.frame >= 0:  # start from arbitrary given reference frame
            COM = aL.center_of_mass()
            COORDs_Li.append(Li.positions - COM)
            COORDs_TFSI.append(TFSI.positions - COM)
    
    dt = np.linspace(0, 5000, 5001)
    msd_Li = finalGetdx_i(COORDs_Li, 0, 1, 24)
    msd_TFSI = finalGetdx_i(COORDs_TFSI, 0, 1, 24)
    msdPrint('msd_m', dt, msd_Li)
    msdPrint('msd_a', dt, msd_TFSI)
    print("MSD analysis finished!")

if __name__ == "__main__":
    main()