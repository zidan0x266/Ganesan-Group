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
import random


def sqPrint(filename, q, sq):  # print sq to the output
    with open('{}.xvg'.format(filename), 'w') as sqout:
        print('# q, S(q)', file=sqout)
        for i in range(len(q)):
            print('{:10.3f} {:10.5f}'.format(q[i], sq[i]), file=sqout)


def qgenerator(cell):
    """
    This function generates the number of q values, since the number
    is depend on the box length (cell), the structure is similar to 
    the following qgenerator function
    """
    counter = 0
    n_point=100000
    n_max = 100
    maxq = 2 * np.pi / 3.0
    minq = 2 * np.pi / cell
    q_vector = np.zeros((n_point, 4))
    for i in range(n_point):
        rand1 = random.random()
        rand2 = random.random()
        rand3 = random.random()
        nx =(n_max + 1) * rand1
        ny =((2 * n_max) + 1) * rand2 - n_max
        nz =((2 * n_max) + 1) * rand3 - n_max
        tmp = (nx * nx + ny * ny + nz * nz)
        if tmp <= (n_max ** 2):
            qx = nx * minq
            qy = ny * minq
            qz = nz * minq
            q2 = qx ** 2 + qy ** 2 + qz ** 2
            q_vector[i] = np.sqrt(q2), qx, qy, qz
    return q_vector


def main():
    uta = mda.Universe("cg_topol.tpr", "cg_traj.xtc")
    # predefined scatter factor for individual atom (atom mass, scatter length)
    scattfactor = dict([(27, 1.0), (28, 1.0), (81, 1.0), 
                        (43, 1.0), (145, 1.0),
                        (41, 0.0), (42, 0.0), (59, 0.0)])
    # initialization
    frequence = 1000
    cell = uta.dimensions[0]  # for cubic box only fetch the x-box length
    q_vector = qgenerator(cell)
    num_q = len(q_vector)
    sq = np.zeros(len(q_vector))
    Counter = np.zeros(len(q_vector))  # counter for normalization
    firstframe = True
    sumscattfact = 0
    print("Entering the S(q) analysis...")
    for frame in uta.trajectory:  # go through all frames
        if frame.frame % frequence == 0:
            print("Currently processing the frame {} ...".format(frame.frame))
            if firstframe:
                factor = np.zeros(frame.n_atoms)
                for aindex in range(frame.n_atoms):  # go through all atoms
                    atomass = int(round(uta.atoms.masses[aindex]))
                    sumscattfact += scattfactor[atomass] ** 2
                    factor[aindex] = scattfactor[atomass]
                firstframe = False
            for q_num in range(len(q_vector)):
                cossum, sinsum = 0, 0
                q_tmp, qx, qy, qz = q_vector[q_num]
                qr = qx * frame.positions[:,0] + qy * frame.positions[:,1] + qz * frame.positions[:,2]
                cossum = np.sum(factor * np.cos(qr))
                sinsum = np.sum(factor * np.sin(qr))
                #for aindex in range(frame.n_atoms):  # go through all atoms
                #    atomass = int(round(uta.atoms.masses[aindex]))
                #    cossum += scattfactor[atomass] * np.cos(qr[aindex])
                #    sinsum += scattfactor[atomass] * np.sin(qr[aindex])
                tmp = np.int(q_tmp / 0.1)
                sq[tmp] += (cossum ** 2 + sinsum ** 2)
                Counter[tmp] += 1
        if frame.frame == 10000:
            break
    q_width = 0.1
    q_vals = []
    Sq_vals = []
    for i in range(1, 400):
        if Counter[i] != 0:
            q_vals.append(np.float(i) * q_width)
            sq[i] /= (sumscattfact * Counter[i])
            Sq_vals.append(sq[i])  # normalization by the number of frames
    print("Output Sq to the data file!")
    sqPrint('sq', q_vals, Sq_vals)


if __name__ == "__main__":
    main()