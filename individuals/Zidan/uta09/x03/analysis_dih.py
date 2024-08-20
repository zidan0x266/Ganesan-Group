"""
Copyright (C) 2018-2022 Zidan Zhang <zhangzidan@gmail.com>
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
from operator import itemgetter


def dihPrint(filename, acf):  # print sqt to the output
    time = np.linspace(0, len(acf) - 1, len(acf))
    with open('{}.xvg'.format(filename), 'w') as sqtout:
        print('# time, dihedral auto-correlation function', file=sqtout)
        for i in range(len(time)):
            print('{:10.3f} {:10.5f}'.format(time[i], acf[i]), file=sqtout)


def getdx_i(DIHEs, t0, t):
    """
    This function give the intermitent Association Lifetime
    C_i = (<cos(t0)cos(t)> - <cos(t0)>**2)/(<cos(t0)cos(t0)> - <cos(t0)>**2) between t0 and t
    """
    avg1 = np.average([i ** 2 for i in DIHEs[t0]])
    avg2 = np.average(DIHEs[t0]) ** 2
    avg3 = np.average([i * j for i, j in zip(DIHEs[t0], DIHEs[t])])
    C_i = (avg3 - avg2) / (avg1 - avg2)
    return C_i


def intervdx_i(DIHEs, t0, tf, reservoir, dt):
    """
    This function gets all the data for the h(t0)h(t0+dt), where
    t0 = 1,2,3,...,tf. This function give us a point of the final plot
    C(t) vs t,
    """
    a = 0
    count = 0
    for t0 in reservoir:
        if t0 + dt <= tf:
            a += getdx_i(DIHEs, t0, t0 + dt)
            count += 1
        else:
            break
    return a / count


def finalGetdx_i(DIHEs, t0, trestart, nt):
    """
    This function gets the final data of the C_i graph.
    """
    tf = len(DIHEs) - 1
    num_points = len(DIHEs)
    pool = mp.Pool(processes=nt)
    reservoir = [x * trestart for x in range(len(DIHEs)) if x * trestart <= tf]
    func = functools.partial(intervdx_i, DIHEs, t0, tf, reservoir)
    return_data = pool.map(func, list(range(num_points)))
    return return_data


def main():
    uta = mda.Universe('../md.tpr', 'md.xtc')
    nframes = len(uta.trajectory)
    backbone = uta.select_atoms('name C01 or name C02')
    dihindex = []
    for i, dih in enumerate(backbone.dihedrals):
        if dih.atoms[0].name == 'C01' and dih.atoms[1].name == 'C02' and dih.atoms[2].name == 'C01' and dih.atoms[-1].name == 'C02':
            dihindex.append(i)
        if dih.atoms[0].name == 'C02' and dih.atoms[1].name == 'C01' and dih.atoms[2].name == 'C02' and dih.atoms[-1].name == 'C01':
            dihindex.append(i)
    cosvals = np.zeros((nframes, len(dihindex)))
    for ts in uta.trajectory:
        frame = ts.frame
        phi = itemgetter(dihindex)(backbone.dihedrals.values())
        cosvals[frame] = np.cos(phi)
        if frame % 10 == 0:
            print('Now processing percentage: {:5.3f}'.format(float(frame) / (nframes - 1)))
    acfphi = finalGetdx_i(cosvals, 0, 10, 56)
    dihPrint('acf_dih', acfphi)


if __name__ == "__main__":
    main()
