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

import argparse
import multiprocessing as mp
import functools

import MDAnalysis as mda
import numpy as np
import logging


logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)-8s %(message)s')
logger = logging.getLogger(__name__)


def _args():
    args = argparse.ArgumentParser("Calculate SQT")
    args.add_argument("tpr")
    args.add_argument("xtc")
    args.add_argument("output")

    return args.parse_args()


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
    largest_q = 1.05
    maximum_q = 40
    num_qvec = np.zeros(maxnumq)
    minq = 2 * np.pi / cell
    for nx in range(-n_max, n_max, 1):
        for ny in range(-n_max, n_max, 1):
            for nz in range(-n_max, n_max, 1):
                if nx == 0 and ny == 0 and nz == 0:
                    continue
                qx = nx * minq; qy = ny * minq; qz = nz * minq
                q_sqrt = np.sqrt(qx**2 + qy**2 + qz**2)
                if q_sqrt > largest_q:
                    continue
                remain = (int(q_sqrt * 1000)) % 100
                if remain >= 90 or remain <= 10:
                    qvalue_10 = int(round(q_sqrt * 10))  # q values multiply by 10
                    if num_qvec[qvalue_10] < maximum_q:                    
                        counter += 1
                        num_qvec[qvalue_10] += 1
    return counter


def qgenerator(n_max, number_q, cell):
    """
    This function generates the q vector
    """
    total_q = number_q
    counter = 0
    largest_q = 1.05
    maximum_q = 40
    num_qvec = np.zeros(total_q)
    q = np.zeros(total_q)
    q_vec = np.zeros((total_q, 3))
    minq = 2 * np.pi / cell  # pi / half_box_length
    for nx in range(-n_max, n_max, 1):
        for ny in range(-n_max, n_max, 1):
            for nz in range(-n_max, n_max, 1):
                if nx == 0 and ny == 0 and nz == 0:
                    continue
                qx = nx * minq; qy = ny * minq; qz = nz * minq
                q_sqrt = np.sqrt(qx**2 + qy**2 + qz**2)
                if q_sqrt > largest_q:
                    continue
                remain = (int(q_sqrt * 1000)) % 100  # assure the generated q has a variance of +- 10%
                if remain >= 90 or remain <= 10:
                    qvalue_10 = int(round(q_sqrt * 10))  # q values multiply by 10
                    if num_qvec[qvalue_10] < maximum_q:                    
                        q_vec[counter][0] = qx
                        q_vec[counter][1] = qy
                        q_vec[counter][2] = qz
                        q[counter] = q_sqrt
                        counter += 1
                        num_qvec[qvalue_10] += 1
    return q_vec, q, counter


def process_frame(qs, uta, scatter_length, frame_idx):
    frame = uta.trajectory[frame_idx]
    poly = uta.select_atoms("resname HPSI or resname MPSI or resname TPSI")
    qr = (qs * poly.positions).sum(axis=1)
    cossum = (scatter_length * np.cos(qr)).sum()
    sinsum = (scatter_length * np.sin(qr)).sum()
    return cossum, sinsum


def main():
    args = _args()

    uta = mda.Universe(args.tpr, args.xtc)
    # only select atoms on the polymer
    poly = uta.select_atoms("resname HPSI or resname MPSI or resname TPSI")
    logger.info(f"Analysing {args.xtc} with topology from {args.xtc}")
    logger.info(f"Frames: {uta.dimensions} {uta.trajectory}")

    # predefined scatter factor for individual atom (atom mass, scatter length)
    # obtained from https://www.ncnr.nist.gov/resources/n-lengths/list.html
    scatter_factor = dict([(1, -3.7406), (12, 6.6511), 
                           (14, 9.3700), (16, 5.8030),
                           (19, 5.6540), (7, -1.9000), 
                           (32, 2.8040)])
    # initialization
    nframes = len(uta.trajectory)
    cell = uta.dimensions[0]  # for cubic box only fetch the x-box length
    number_q = preq(21, cell)
    q_vec, q, n_max = qgenerator(21, number_q, cell)
    sqt = np.zeros((n_max, nframes))
    sqtvals = np.zeros((11, nframes))  # bins from 0.0 to 1.0 with increament 0.1
    sqttmps = np.zeros((11, nframes))

    scatter_length = np.array(list(map(scatter_factor.get, poly.masses.round().astype(int))))

    pool = mp.Pool(20)

    for qindex in range(n_max):  # go through all q values
        cossum = np.zeros(nframes)
        sinsum = np.zeros(nframes)
        frame_indexes = np.arange(nframes, dtype = int)
        fn_frame_analysis = functools.partial(process_frame, q_vec[qindex], uta, scatter_length)
        logger.info(f"Start processing frames for {qindex}")
        frame_map = pool.map(fn_frame_analysis, frame_indexes)

        for tidx, f_sum in enumerate(frame_map):  # get returned fourier transform summation
            cossum[tidx] = f_sum[0]
            sinsum[tidx] = f_sum[1]
        logger.info(f"Finished processing frame")
        for i in range(nframes):
            for j in range(0, nframes - i):
                k = i + j
                sqt[qindex][i] += (2.0 * (cossum[j] * cossum[k] + sinsum[j] * sinsum[k]))
        # averaging for the statistic enhancement
        for i in range(nframes):
            sqt[qindex][i] /= (2.0 * (nframes - i))
    logger.info("The raw data are prepared successfully!")
    # preparing the final sqt values with bining average
    logger.info("Entering the bin averaging!")
    for i in range(len(sqttmps)):
        for qindex in range(n_max):
            if round(q[qindex] * 10) == i:
                for time in range(nframes):
                    sqttmps[i][time] += sqt[qindex][time]
    # normallization
    logger.info("Entering the normalization!")
    for i in range(len(sqtvals)):
        if sqttmps[i][0] != 0:
            sqtvals[i] = sqttmps[i] / sqttmps[i][0]
    logger.info("Output Sq(t) to the data file!")
    sqtPrint(args.output, sqtvals)


if __name__ == "__main__":
    main()
