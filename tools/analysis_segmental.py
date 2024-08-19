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


def sqtPrint(filename, sqt):  # print sqt to the output
    with open('{}.xvg'.format(filename), 'w') as sqtout:
        print('# time 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0', file=sqtout)
        for i in range(len(sqt[0])):
            print('{:10.3f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(i, 
            sqt[0][i], 
            sqt[1][i], 
            sqt[2][i], 
            sqt[3][i], 
            sqt[4][i], 
            sqt[5][i], 
            sqt[6][i], 
            sqt[7][i], 
            sqt[8][i], 
            sqt[9][i], 
            sqt[10][i]), file=sqtout)


def preq(n_max, cell):
    """
    This function generates the number of q values, since the number
    is depend on the box length (cell), the structure is similar to 
    the following qgenerator function
    """
    maxnumq = 1000
    counter = 0
    N = np.zeros(maxnumq)
    minq = 2 * np.pi / cell
    for nx in range(-n_max, n_max, 1):
        for ny in range(-n_max, n_max, 1):
            for nz in range(-n_max, n_max, 1):
                if nx == 0 and ny == 0 and nz == 0:
                    continue
                qx = nx * minq; qy = ny * minq; qz = nz * minq
                q = np.sqrt(qx**2 + qy**2 + qz**2)
                if q > 1.05:
                    continue
                remain = (int(q * 1000)) % 100
                if remain >= 90 or remain <= 10:
                    a = int(round(q * 10))
                    if N[a] < 40:                    
                        counter += 1
                        N[a] += 1
    return counter


def qgenerator(n_max, counter, cell):
    """
    This function generates the q vector
    """
    maxnumq = counter
    counter = 0
    N = np.zeros(maxnumq)
    N2 = np.zeros(maxnumq)
    q = np.zeros(maxnumq)
    n_vec = np.zeros((maxnumq, 3))
    minq = 2 * np.pi / cell
    for nx in range(-n_max, n_max, 1):
        for ny in range(-n_max, n_max, 1):
            for nz in range(-n_max, n_max, 1):
                if nx == 0 and ny == 0 and nz == 0:
                    continue
                qx = nx * minq; qy = ny * minq; qz = nz * minq
                q2 = np.sqrt(qx**2 + qy**2 + qz**2)
                if q2 > 1.05:
                    continue
                remain = (int(q2 * 1000)) % 100
                if remain >= 90 or remain <= 10:
                    a = int(round(q2 * 10))
                    if N2[a] < 40:                    
                        n_vec[counter][0] = qx
                        n_vec[counter][1] = qy
                        n_vec[counter][2] = qz
                        q[counter] = q2
                        counter += 1
                        N2[a] += 1
    #print(counter)
    return n_vec, q, counter


def main():
    uta = mda.Universe("md.tpr", "unwrap.xtc")
    # predefined scatter factor for individual atom (atom mass, scatter length)
    scattfactor = dict([(1, -3.7406), (12, 6.6511), 
                        (14, 9.3700), (16, 5.8030),
                        (19, 5.6540), (7, -1.9000), 
                        (32, 2.8040)])
    # initialization
    nframes = len(uta.trajectory)
    cell = uta.dimensions[0]  # for cubic box only fetch the x-box length
    numberq = preq(21, cell)
    n_vec, q, n_max = qgenerator(21, numberq, cell)
    sqt = np.zeros((n_max, nframes))
    sqtvals = np.zeros((11, nframes))  # bins from 0.0 to 1.0 with increament 0.1
    sqttmps = np.zeros((11, nframes))
    # select heavy atoms on the polymer
    aS = uta.select_atoms('name C1 or name C2 or name N1 or name N2 or name C3 or name C4 or name C5 or name C6 or name C7 or name C8 or name C9')
    for qindex in range(n_max):  # go through all q values
        cossum = np.zeros(nframes)
        sinsum = np.zeros(nframes)
        qx = n_vec[qindex][0]
        qy = n_vec[qindex][1]
        qz = n_vec[qindex][2]
        for frame in uta.trajectory:  # go through all frames
            for aindex in range(len(aS)):  # go through all atoms
                qr = qx * aS.positions[aindex][0] + qy * aS.positions[aindex][1] + qz * aS.positions[aindex][2]
                atomass = int(round(aS.masses[aindex]))
                time = frame.frame
                cossum[time] += scattfactor[atomass] * np.cos(qr)
                sinsum[time] += scattfactor[atomass] * np.sin(qr)
        for i in range(nframes):
            for j in range(0, nframes - i):
                k = i + j
                sqt[qindex][i] += (2.0 * (cossum[j] * cossum[k] + sinsum[j] * sinsum[k]))
        # averaging for the statistic enhancement
        for i in range(nframes):
            sqt[qindex][i] /= (2.0 * (nframes - i))
    print("The raw data are prepared successfully!")
    # preparing the final sqt values with bining average
    print("Entering the bin averaging!")
    for i in range(len(sqttmps)):
        for qindex in range(n_max):
            if round(q[qindex] * 10) == i:
                for time in range(nframes):
                    sqttmps[i][time] += sqt[qindex][time]
    # normallization
    print("Entering the normalization!")
    for i in range(len(sqtvals)):
        if sqttmps[i][0] != 0:
            sqtvals[i] = sqttmps[i] / sqttmps[i][0]
    print("Output Sq(t) to the data file!")
    sqtPrint('sqt', sqtvals)


if __name__ == "__main__":
    main()