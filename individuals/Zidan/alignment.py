#!/usr/bin/env python
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

vex = [1, 0, 0]
vey = [0, 1, 0]
vez = [0, 0, 1]

uta = mda.Universe('topol.tpr', 'unwrap.xtc')
BB = uta.select_atoms('type BB')
en_x = []; en_y = []; en_z = []
for ts in uta.trajectory:
    chain = []
    for i in range(200):
        single = []
        for j in range(8):
            k = i * 8 + j
            single.append(BB.positions[k])
        chain.append(single)

    chain_vec = []
    for i in range(200):
        tmp = []
        for j in range(7):
            tmp.append(chain[i][j+1] - chain[i][j])
        chain_vec.append(chain[i][-1] - chain[i][0])
    
    for i in range(200):
        en_x.append(angle_between(vex, chain_vec[i]))
        en_y.append(angle_between(vey, chain_vec[i]))
        en_z.append(angle_between(vez, chain_vec[i]))

xx = np.cos(en_x)**2
yy = np.cos(en_y)**2
zz = np.cos(en_z)**2

s_x = [0.5 * (3 * xx[i] - 1) for i in range(len(xx))]
s_y = [0.5 * (3 * yy[i] - 1) for i in range(len(yy))]
s_z = [0.5 * (3 * zz[i] - 1) for i in range(len(zz))]
#print(np.average(en_x), np.std(en_x))
#print(np.average(en_y), np.std(en_y))
#print(np.average(en_z), np.std(en_z))
print(np.average(s_x), np.std(s_x))
print(np.average(s_y), np.std(s_y))
print(np.average(s_z), np.std(s_z))
