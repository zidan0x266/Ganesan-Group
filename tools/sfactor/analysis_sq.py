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
import random
import logging


logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)-8s %(message)s')
logger = logging.getLogger(__name__)


def _args():
    args = argparse.ArgumentParser("Calculate Static SQ")
    args.add_argument("tpr")
    args.add_argument("xtc")
    args.add_argument("output")

    return args.parse_args()


def sqPrint(filename, sq):  # print sq to the output
    with open('{}.xvg'.format(filename), 'w') as sqtout:
        print('# q S(q)', file=sqtout)
        for i in range(len(sq[0])):
            print('{:10.3f} {:10.5f}'.format(sq[0][i], sq[1][i]), file=sqtout)


def qgenerator(number_q, cell):
    """
    This function generates the q vector
    """
    total_q = number_q
    counter = 0
    RAND_MAX = 2 ** 15 - 1  # from C++
    n_max = 100
    num_qvec = np.zeros(total_q)
    q = np.zeros(total_q)
    q_vec = np.zeros((total_q, 3))
    minq = 2 * np.pi / cell  # pi / half_box_length
    for i in range(total_q):
        newborn = True
        while newborn:
            rand_x = random.randrange(RAND_MAX + 1) / RAND_MAX
            rand_y = random.randrange(RAND_MAX + 1) / RAND_MAX
            rand_z = random.randrange(RAND_MAX + 1) / RAND_MAX
            nx =(n_max + 1) * rand_x
            ny =((2 * n_max) + 1) * rand_y - n_max
            nz =((2 * n_max) + 1) * rand_z - n_max
            n2 = nx ** 2 + ny ** 2 + nz ** 2
            if n2 <= n_max ** 2:
                qx = nx * minq
                qy = ny * minq
                qz = nz * minq
                newborn = False
        q2 = qx**2 + qy**2 + qz**2
        q[i] = np.sqrt(q2)
        q_vec[i][0] = qx
        q_vec[i][1] = qy
        q_vec[i][2] = qz
    return q_vec, q


def process_frame(uta, atom_group, cell, number_q, q_width, scatter_length, frame_idx):
    frame = uta.trajectory[frame_idx]
    poly = uta.select_atoms(atom_group)
    q_vec, q = qgenerator(number_q, cell)
    q_points = int(1 / q_width) * 10 + 1
    sq_frame = np.zeros(q_points)
    N_frame = np.zeros(q_points)
    for qindex in range(number_q):
        qr = (q_vec[qindex] * poly.positions).sum(axis=1)
        cossum = (scatter_length * np.cos(qr)).sum()
        sinsum = (scatter_length * np.sin(qr)).sum()
        qvalue_10 = int(q[qindex] * int(1 / q_width))
        sq_frame[qvalue_10] += cossum ** 2 + sinsum ** 2
        N_frame[qvalue_10] += 1
    return sq_frame, N_frame


def main():
    args = _args()

    uta = mda.Universe(args.tpr, args.xtc)
    # only select atoms on the polymer
    atom_group = "resname HPSI or resname MPSI or resname TPSI"
    poly = uta.select_atoms(atom_group)
    logger.info(f"Analysing {args.xtc} with topology from {args.tpr}")
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
    number_q = 10000
    q_width = 0.1
    q_points = int(1 / q_width) * 10 + 1
    N = np.zeros((nframes, q_points))
    sqtmps = np.zeros((nframes, q_points))
    sqvals = np.zeros((2, q_points))  # bins from 0.0 to 1.0 with increament 0.1

    scatter_length = np.array(list(map(scatter_factor.get, poly.masses.round().astype(int))))
    # scatter_length_sum corresponds to the summation of bi*bj for the normalization
    # eq.2 in Ref: https://pubs.acs.org/doi/full/10.1021/acs.macromol.7b00724
    scatter_length_sum = np.square(scatter_length).sum()

    pool = mp.Pool(20)
    frame_indexes = np.arange(nframes, dtype = int)
    frame_analysis = functools.partial(process_frame, uta, atom_group, cell, number_q, q_width, scatter_length)
    frame_map = pool.map(frame_analysis, frame_indexes)
    for tidx, sq_per_frame in enumerate(frame_map):  # get returned fourier transform summation
        sqtmps[tidx] = sq_per_frame[0]
        N[tidx] = sq_per_frame[1]
    logger.info(f"Finished processing frame")

    logger.info("Entering the normalization!")
    sq_sum = np.sum(sqtmps, axis = 0)
    N_sum = np.sum(N, axis = 0)
    for i in range(1, q_points):
        sqvals[0][i] = i * q_width
        if N_sum[i] != 0:
            sqvals[1][i] = sq_sum[i] / (scatter_length_sum * N_sum[i])
    logger.info("Output S(q) to the data file!")
    sqPrint(args.output, sqvals)


if __name__ == "__main__":
    main()
