"""
Copyright (C) 2018-2023 Zidan Zhang <zhangzidan@gmail.com>

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

import functools
import argparse
import numpy as np
import random
import logging


logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)-8s %(message)s')
logger = logging.getLogger(__name__)


class MyArgParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(MyArgParser, self).__init__(*args, **kwargs)

    def convert_arg_line_to_args(self, line):
        for arg in line.split():
            t = arg.strip()
            if not t:
                continue
            if t.startswith('#'):
                break
            if not t.startswith('--'):
                t = '--{}'.format(t)
            yield t

    @staticmethod
    def save_to_file(output_file, namespace):
        """Saves arguments to file so it can be read again.

        Args:
            output_file: The string with the name of output file.
            namespace: The namespace with arguments.
        """
        with open(output_file, "w") as of:
            for k, v in namespace.__dict__.iteritems():
                if v is not None:
                    of.write('{}={}\n'.format(k, v))


def _args():
    parser = MyArgParser(description="Calculate Static SQ", fromfile_prefix_chars='@')
    parser.add_argument("--xyz", type=str, help="input coordinates file with atomtype information")
    parser.add_argument("--xvg", default="sq.xvg", type=str, help="output S(q) to .xvg file")

    return parser.parse_args()


def sqPrint(filename, sq):  # print sq to the output
    with open('{}'.format(filename), 'w') as sqtout:
        print('# q S(q)', file=sqtout)
        for i in range(len(sq[0])):
            print('{:10.3f} {:10.5f}'.format(sq[0][i], sq[1][i]), file=sqtout)


def qgenerator(number_q, cell):
    """
    This function generates the q vector
    """
    total_q = number_q
    RAND_MAX = 2 ** 15 - 1  # from C++
    n_max = 100
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


def main():
    args = _args()

    uta = np.loadtxt(args.xyz)
    # only select atoms on the polymer
    atoms = [int(x) for x in uta[:,0]]
    positions = uta[:,1:]
    logger.info(f"Analysing {args.xyz} for the reactive Monte-Carlo")

    # predefined scatter factor for individual atom (atom mass, scatter length)
    # obtained from https://www.ncnr.nist.gov/resources/n-lengths/list.html
    scatter_factor = dict([(1, -3.7406), (2, 6.6511), (3, 2.8040)])
    # initialization
    cell = 50  # for cubic box only fetch the x-box length
    number_q = 20000
    q_width = 0.1
    q_points = int(1 / q_width) * 10 + 1
    N = np.zeros(number_q)
    sqtmps = np.zeros(number_q)
    sqvals = np.zeros((2, q_points))  # bins from 0.0 to 1.0 with increament 0.1

    scatter_length = np.array(list(map(scatter_factor.get, np.array(atoms).astype(int))))
    # scatter_length_sum corresponds to the summation of bi*bj for the normalization
    # eq.2 in Ref: https://pubs.acs.org/doi/full/10.1021/acs.macromol.7b00724
    scatter_length_sum = np.square(scatter_length).sum()

    q_vec, q = qgenerator(number_q, cell)
    for qindex in range(number_q):
        qr = (q_vec[qindex] * positions).sum(axis=1)
        cossum = (scatter_length * np.cos(qr)).sum()
        sinsum = (scatter_length * np.sin(qr)).sum()
        qvalue_10 = int(q[qindex] * 10)
        sqtmps[qvalue_10] += cossum ** 2 + sinsum ** 2
        N[qvalue_10] += 1
    logger.info("Entering the normalization!")
    for i in range(q_points):
        sqvals[0][i] = i * q_width
        if N[i] != 0:
            sqvals[1][i] = sqtmps[i] / (scatter_length_sum * N[i])
    logger.info("Output S(q) to the data file!")
    sqPrint(args.xvg, sqvals)


if __name__ == "__main__":
    main()
