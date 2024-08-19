#!/usr/bin/env python3
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

import MDAnalysis as mda
import numpy as np

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    if angle / np.pi <= 0.5:
        return angle
    else:
        return np.pi - angle


def main():
    vex = [1, 0, 0]
    vey = [0, 1, 0]
    vez = [0, 0, 1]
    top = 'cg_topol.tpr'
    trj = 'cgtrj.xtc'
    Nchains = 100
    chainlength = 20
    uta = mda.Universe(top, trj)
    BB = uta.select_atoms('type BB')
    en_x = []; en_y = []; en_z = []
    for ts in uta.trajectory:
        chain = []
        for i in range(Nchains):
            single = []
            for j in range(chainlength):
                k = i * chainlength + j
                single.append(BB.positions[k])
            chain.append(single)    
        chain_vec = []
        sd1 = []; sd2 = []; sd3 = []
        for i in range(Nchains):
            tmp = []; tmp_x = []; tmp_y = []; tmp_z = []
            for j in range(chainlength - 1):
                tmp = chain[i][j+1] - chain[i][j]
                tmp_x.append(angle_between(vex, tmp))
                tmp_y.append(angle_between(vey, tmp))
                tmp_z.append(angle_between(vez, tmp))
            xx = np.cos(tmp_x)**2
            yy = np.cos(tmp_y)**2
            zz = np.cos(tmp_z)**2
            s_x = [0.5 * (3 * xx[i] - 1) for i in range(len(xx))]
            s_y = [0.5 * (3 * yy[i] - 1) for i in range(len(yy))]
            s_z = [0.5 * (3 * zz[i] - 1) for i in range(len(zz))]
            sd1.append(np.average(s_x))
            sd2.append(np.average(s_y))
            sd3.append(np.average(s_z))           
        en_x.append(np.average(sd1))
        en_y.append(np.average(sd2))
        en_z.append(np.average(sd3))    
        if ts.frame % 10 == 0:
            print("Frame: {0:5d}, Time: {1:8.3f} ps".format(ts.frame, uta.trajectory.time))
        if ts.frame == 100:
            break   
    print(np.average(en_x), np.std(en_x))
    print(np.average(en_y), np.std(en_y))
    print(np.average(en_z), np.std(en_z))


if __name__ == "__main__":
    main()