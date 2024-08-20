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

import MDAnalysis as mda
import numpy as np


def dpPrint(filename, dt, dx):  # print displacement
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# time delta-X', file=anaout)
        for i in range(0, len(dt)):
            print('{:10.3f} {:10.5f}'.format(dt[i], dx[i]), file=anaout)


def getdx_i(COORDs, t0, t):
    C_i = COORDs[t] - COORDs[t0]
    return np.sum(np.square(C_i), axis = 1)


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
    tf = len(COORDs) - 1
    num_points = len(COORDs)
    pool = mp.Pool(processes=nt)
    reservoir = [x * trestart for x in range(tf) if x * trestart < tf]
    func = functools.partial(intervdx_i, COORDs, tf, reservoir)
    return_data = pool.map(func, list(range(num_points)))
    return return_data


def main():
    top = "../cg_topol.tpr"
    trj = "unwrap.xtc"
    trestart = 10
    nproc = 80
    withsalt = False
    stratframe = 0
    uta = mda.Universe(top, trj)
    t_series = np.linspace(0, len(uta.trajectory) - 1, len(uta.trajectory))
    if withsalt:
        gL = uta.select_atoms('all')
        g1 = uta.select_atoms('type TF')
        g2 = uta.select_atoms('type IM')
        g3 = uta.select_atoms('type LI')
        g1_coords, g2_coords, g3_coords = [], [], []
        for ts in uta.trajectory:
            if ts.frame >= stratframe:
                ref = gL.center_of_mass()
                g1_coords.append(g1.positions - ref)
                g2_coords.append(g2.positions - ref)
                g3_coords.append(g3.positions - ref)
                if ts.frame % 100 == 0:
                    print("Frame {}".format(ts.frame))
        dxg1 = finalGetdx_i(g1_coords, trestart, nproc)
        disp_g1 = np.average(dxg1, axis = 1)
        dpPrint('msd_TF', t_series, disp_g1)
        dxg2 = finalGetdx_i(g2_coords, trestart, nproc)
        disp_g2 = np.average(dxg2, axis = 1)
        dpPrint('msd_IM', t_series, disp_g2)
        dxg3 = finalGetdx_i(g3_coords, trestart, nproc)
        disp_g3 = np.average(dxg3, axis = 1)
        dpPrint('msd_LI', t_series, disp_g3)
        print("Displacement analysis finished!")
    else:
        #gL = uta.select_atoms('all')  # System center of mass
        gL = uta.select_atoms('type BB or type BZ or type IM')  # Polymer chain center of mass
        g1 = uta.select_atoms('type TF')
        g2 = uta.select_atoms('type IM')
        g1_coords, g2_coords = [], []
        for ts in uta.trajectory:
            if ts.frame >= stratframe:
                ref = gL.center_of_mass()
                g1_coords.append(g1.positions - ref)
                g2_coords.append(g2.positions - ref)
                if ts.frame % 100 == 0:
                    print("Frame {}".format(ts.frame))
        dxg1 = finalGetdx_i(g1_coords, trestart, nproc)
        disp_g1 = np.average(dxg1, axis = 1)
        dpPrint('msd_TF_POM', t_series, disp_g1)
        dxg2 = finalGetdx_i(g2_coords, trestart, nproc)
        disp_g2 = np.average(dxg2, axis = 1)
        dpPrint('msd_IM_POM', t_series, disp_g2)
        print("MSD analysis finished!")
    

if __name__ == "__main__":
    main()