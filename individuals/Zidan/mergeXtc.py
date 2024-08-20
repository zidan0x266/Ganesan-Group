"""
Copyright (C) 2018-2021 Zidan Zhang <zhangzidan@gmail.com>
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

import MDAnalysis
import MDAnalysis as mda


def main():
    top = "cg_topol.tpr"
    trj1 = "cg_traj_01.xtc"
    trj2 = "cg_traj_02.xtc"
    uta = mda.Universe(top, trj1, trj2, continuous=True)
    aL = uta.select_atoms("all")
    with MDAnalysis.Writer("data.xtc", aL.n_atoms) as W:
        for ts in uta.trajectory:
            if ts.frame != 10000:
                W.write(aL)


if __name__ == "__main__":
    main()