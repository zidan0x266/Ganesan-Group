"""
Copyright (C) 2018-2022 Zidan Zhang <zhangzidan@gmail.com>
​
This file is distributed under free software licence:
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
​
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
​
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import seaborn as sn
import os
from os import listdir
from os.path import isfile, join
from scipy import stats


def loadMSD(filepath, length, func, Natoms):
    if func == 5:
        xdata = True
        RAW = np.zeros((length, func + 1))
        MSD = np.zeros((length, func + 1))
        rawfiles = [f for f in listdir(filepath) if isfile(join(filepath, f))]
        for f in rawfiles:
            if f.startswith('msd_self_TF'):
                temp = np.loadtxt(filepath + f)
                if xdata:
                    RAW[:, 0] = temp[:, 0]
                    xdata = False
                RAW[:, 1] = temp[:, 1]
            elif f.startswith('msd_self_IM'):
                temp = np.loadtxt(filepath + f)
                if xdata:
                    RAW[:, 0] = temp[:, 0]
                    xdata = False
                RAW[:, 2] = temp[:, 1]
            elif f.startswith('msd_TF_TF'):
                temp = np.loadtxt(filepath + f)
                if xdata:
                    RAW[:, 0] = temp[:, 0]
                    xdata = False
                RAW[:, 3] = temp[:, 1]
            elif f.startswith('msd_IM_IM'):
                temp = np.loadtxt(filepath + f)
                if xdata:
                    RAW[:, 0] = temp[:, 0]
                    xdata = False
                RAW[:, 4] = temp[:, 1]
            elif f.startswith('msd_TF_IM'):
                temp = np.loadtxt(filepath + f)
                if xdata:
                    RAW[:, 0] = temp[:, 0]
                    xdata = False
                RAW[:, 5] = temp[:, 1]
        MSD[:, 0] = RAW[:, 0] / 1000  # Convert ps to ns
        MSD[:, 1] = RAW[:, 1] / Natoms[0] / 100
        MSD[:, 2] = RAW[:, 2] / Natoms[1] / 100
        MSD[:, 3] = (RAW[:, 3] - RAW[:, 1]) / Natoms[0] / 100
        MSD[:, 4] = (RAW[:, 4] - RAW[:, 2]) / Natoms[1] / 100
        MSD[:, 5] = RAW[:, 5] / ((Natoms[0] * Natoms[1]) / (Natoms[0] + Natoms[1])) / 100
    elif func == 9:
        xdata = True
        RAW = np.zeros((length, func + 1))
        MSD = np.zeros((length, func + 1))
        rawfiles = [f for f in listdir(filepath) if isfile(join(filepath, f))]
        for f in rawfiles:
            if f.startswith('msd_self_TF'):
                temp = np.loadtxt(filepath + f)
                if xdata:
                    RAW[:, 0] = temp[:, 0]
                    xdata = False
                RAW[:, 1] = temp[:, 1]
            elif f.startswith('msd_self_IM'):
                temp = np.loadtxt(filepath + f)
                if xdata:
                    RAW[:, 0] = temp[:, 0]
                    xdata = False
                RAW[:, 2] = temp[:, 1]
            elif f.startswith('msd_self_LI'):
                temp = np.loadtxt(filepath + f)
                if xdata:
                    RAW[:, 0] = temp[:, 0]
                    xdata = False
                RAW[:, 3] = temp[:, 1]
            elif f.startswith('msd_TF_TF'):
                temp = np.loadtxt(filepath + f)
                if xdata:
                    RAW[:, 0] = temp[:, 0]
                    xdata = False
                RAW[:, 4] = temp[:, 1]
            elif f.startswith('msd_IM_IM'):
                temp = np.loadtxt(filepath + f)
                if xdata:
                    RAW[:, 0] = temp[:, 0]
                    xdata = False
                RAW[:, 5] = temp[:, 1]
            elif f.startswith('msd_LI_LI'):
                temp = np.loadtxt(filepath + f)
                if xdata:
                    RAW[:, 0] = temp[:, 0]
                    xdata = False
                RAW[:, 6] = temp[:, 1]
            elif f.startswith('msd_TF_IM'):
                temp = np.loadtxt(filepath + f)
                if xdata:
                    RAW[:, 0] = temp[:, 0]
                    xdata = False
                RAW[:, 7] = temp[:, 1]
            elif f.startswith('msd_TF_LI'):
                temp = np.loadtxt(filepath + f)
                if xdata:
                    RAW[:, 0] = temp[:, 0]
                    xdata = False
                RAW[:, 8] = temp[:, 1]
            elif f.startswith('msd_IM_LI'):
                temp = np.loadtxt(filepath + f)
                if xdata:
                    RAW[:, 0] = temp[:, 0]
                    xdata = False
                RAW[:, 9] = temp[:, 1]
        MSD[:, 0] = RAW[:, 0] / 1000  # Convert ps to ns
        MSD[:, 1] = RAW[:, 1] / Natoms[0] / 100
        MSD[:, 2] = RAW[:, 2] / Natoms[1] / 100
        MSD[:, 3] = RAW[:, 3] / Natoms[2] / 100
        MSD[:, 4] = (RAW[:, 4] - RAW[:, 1]) / Natoms[0] / 100
        MSD[:, 5] = (RAW[:, 5] - RAW[:, 2]) / Natoms[1] / 100
        MSD[:, 6] = (RAW[:, 6] - RAW[:, 3]) / Natoms[2] / 100
        MSD[:, 7] = RAW[:, 7] / ((Natoms[0] * Natoms[1]) / (Natoms[0] + Natoms[1])) / 100
        MSD[:, 8] = RAW[:, 8] / ((Natoms[0] * Natoms[2]) / (Natoms[0] + Natoms[2])) / 100
        MSD[:, 9] = RAW[:, 9] / ((Natoms[1] * Natoms[2]) / (Natoms[1] + Natoms[2])) / 100
    return MSD


def getDiff(MSD, fitrange):
    Diff = []
    for i in range(1, MSD.shape[1]):
        Diff.append(np.average(MSD[fitrange[0]:fitrange[1], i] / MSD[fitrange[0]:fitrange[1], 0]) * 1e-5 / 6)
    return Diff


def main(folder):
    Diff = []
    for sample in ['s01', 's02', 's03', 's04', 's05']:
        path = sample + '/' + folder + '/'
        length = 100001
        func = 5
        Natoms = [2000, 2000, 0]
        fitrange = [20000, 80000]
        data = loadMSD(path, length, func, Natoms)
        temp = getDiff(data, fitrange)
        Diff.append(temp)
    print(np.average(Diff, axis = 0))
    print(np.std(Diff, axis = 0))


if __name__ == "__main__":
    main('c000')