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

"""
Create density grid file using VMD in command line
mol new densmap.gro
volmap density [atomselect top "resname HMMA MMMA TMMA"] -res 0.5 -weight mass -o dA.dx
volmap density [atomselect top "resname HVIM MVIM TVIM PF6"] -res 0.5 -weight mass -o dB.dx
quit
"""

import numpy as np
from gridData import Grid

def outPrint(filename, vol1, vol0, flac):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# vol1, vol0, flac', file=anaout)
        print('{:20.1f} {:20.1f} {:10.5f}'.format(vol1, vol0, flac), file=anaout)


def main():
    dA = Grid("dA.dx")
    dB = Grid("dB.dx")
    dim = np.max((np.max(dA.grid.shape), np.max(dA.grid.shape)))
    print(dim)
    dA_new = np.zeros((dim, dim, dim))
    dB_new = np.zeros((dim, dim, dim))
    for i in range(dA.grid.shape[0]):
        for j in range(dA.grid.shape[1]):
            for k in range(dA.grid.shape[2]):
                dA_new[i][j][k] = dA.grid[i][j][k]
    for i in range(dB.grid.shape[0]):
        for j in range(dB.grid.shape[1]):
            for k in range(dB.grid.shape[2]):
                dB_new[i][j][k] = dB.grid[i][j][k]
    dI = dA_new * dB_new
    out = Grid("dA.dx")
    for i in range(out.grid.shape[0]):
        for j in range(out.grid.shape[1]):
            for k in range(out.grid.shape[2]):
                out.grid[i][j][k] = dI[i][j][k]
    out.export("if.dx")
    vol1 = len(out.grid[np.where(out.grid > 0)]) * 1.0
    vol0 = out.grid.shape[0] * out.grid.shape[1] * out.grid.shape[2] * 1.0
    frac = vol1 / vol0
    outPrint('interfacialV', vol1, vol0, frac)


if __name__ == "__main__":
    main()