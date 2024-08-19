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

import numpy as np
import scipy.linalg as la


def outPrint(filename, eigens):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# lambda_1, lambda_2, lambda_3', file=anaout)
        for i in range(len(eigens)):
            print('{:10.5f} {:10.5f} {:10.5f}'.format(eigens[i][0], eigens[i][1], eigens[i][2]), file=anaout)


def getTensor(msdt):
    tensor = np.zeros((3, 3))
    tensor[0][0] = msdt[2]
    tensor[0][1] = msdt[5]
    tensor[0][2] = msdt[6]
    tensor[1][0] = msdt[5]
    tensor[1][1] = msdt[3]
    tensor[1][2] = msdt[7]
    tensor[2][0] = msdt[6]
    tensor[2][1] = msdt[2]
    tensor[2][2] = msdt[7]
    return tensor


def main():
    MSD = np.loadtxt('msd_PF6.xvg')
    Eigens = []
    for frame in range(MSD.shape[0]):
        msdmatrix = getTensor(MSD[frame])
        eigvals, _ = la.eig(msdmatrix)
        Eigens.append(eigvals.real)
    outPrint('msd_eigens', Eigens)


if __name__ == "__main__":
    main()