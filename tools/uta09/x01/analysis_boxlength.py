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

from os import lseek
import MDAnalysis as mda
import numpy as np

def main():
    top = "gentpr/cg_topol.tpr"
    trj = "cg_traj.xtc"
    stratframe = 20000
    uta = mda.Universe(top, trj)
    length = []
    for ts in uta.trajectory:
        if ts.frame >= stratframe:
            cell = ts.dimensions
            length.append(cell[0])
        if ts.frame % 100 == 0:
            print("Frame {}".format(ts.frame))
    print(np.average(length))


if __name__ == "__main__":
    main()