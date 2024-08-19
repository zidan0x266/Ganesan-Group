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

import glob
import numpy as np

def dpPrint(filename, dt, dx):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# time dx', file=anaout)
        for i in range(0, len(dt)):
            print('{:10.3f} {:10.5f}'.format(dt[i], dx[i]), file=anaout)

filenames = sorted(glob.glob('*xvg'))
samples = []
for f in filenames:
    samples.append(np.loadtxt(f, usecols = 1))
avg = np.average(samples, axis = 0)
dt = np.linspace(0, len(avg) - 1, len(avg))
dpPrint('anion', dt, avg)