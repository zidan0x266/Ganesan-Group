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

import numpy as np
import scipy.linalg as la
import os
from os import listdir
from os.path import isfile, join
from scipy import stats


def loadrawMSD(filepath, fn):
    MSD = []
    rawfiles = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    for f in rawfiles:
        if f.startswith(fn):
            MSD.append(np.loadtxt(filepath + f))
    return MSD


def getDiff_ho(MSD, fitrange):
    Diff = np.average(MSD[fitrange[0]:fitrange[1], 1] / (MSD[fitrange[0]:fitrange[1], 0] / 1000)) * 1e-5 / 6
    return Diff


def getDtensor(MSD, fitrange):
    tensorD = np.zeros((3, 3))
    dxx = np.average(MSD[fitrange[0]:fitrange[1], 2] / (MSD[fitrange[0]:fitrange[1], 0] / 1000)) * 1e-5 / 2
    dyy = np.average(MSD[fitrange[0]:fitrange[1], 3] / (MSD[fitrange[0]:fitrange[1], 0] / 1000)) * 1e-5 / 2
    dzz = np.average(MSD[fitrange[0]:fitrange[1], 4] / (MSD[fitrange[0]:fitrange[1], 0] / 1000)) * 1e-5 / 2
    dyx = np.average(MSD[fitrange[0]:fitrange[1], 5] / (MSD[fitrange[0]:fitrange[1], 0] / 1000)) * 1e-5 / 2
    dzx = np.average(MSD[fitrange[0]:fitrange[1], 6] / (MSD[fitrange[0]:fitrange[1], 0] / 1000)) * 1e-5 / 2
    dzy = np.average(MSD[fitrange[0]:fitrange[1], 7] / (MSD[fitrange[0]:fitrange[1], 0] / 1000)) * 1e-5 / 2
    tensorD[0][0] = dxx
    tensorD[0][1] = dyx
    tensorD[0][2] = dzx
    tensorD[1][0] = dyx
    tensorD[1][1] = dyy
    tensorD[1][2] = dzy
    tensorD[2][0] = dzx
    tensorD[2][1] = dzy
    tensorD[2][2] = dzz
    return tensorD


def main():
    filepath = './'
    file_prefix = ['msd_tf_', 'msd_im_']
    func_type = {0: 'HO', 1: 'LA'}
    func_key = 1
    fitrange = [30000, 39000]
    if func_type[func_key] == 'HO':
        shortDa = []
        shortDc = []
        rawTF = loadrawMSD(filepath, file_prefix[0])
        rawIM = loadrawMSD(filepath, file_prefix[1])
        for i in range(len(rawTF)):
            shortDa.append(getDiff_ho(rawTF[i], fitrange))
            shortDc.append(getDiff_ho(rawIM[i], fitrange))
        print("{:8.5e}".format(np.average(shortDa)))
        print("{:8.5e}".format(np.average(shortDc)))
    elif func_type[func_key] == 'LA':
        shortDa = []
        shortDc = []
        rawTF = loadrawMSD(filepath, file_prefix[0])
        rawIM = loadrawMSD(filepath, file_prefix[1])
        for i in range(len(rawTF)):
            tmp_vals_a, _ = la.eig(getDtensor(rawTF[i], fitrange))
            shortDa.append(np.average([x for x in tmp_vals_a.real if x != np.min(tmp_vals_a.real)]))
            tmp_vals_c, _ = la.eig(getDtensor(rawIM[i], fitrange))
            shortDc.append(np.average([x for x in tmp_vals_c.real if x != np.min(tmp_vals_c.real)]))
        print("{:8.5e}".format(np.average(shortDa)))
        print("{:8.5e}".format(np.average(shortDc)))


if __name__ == "__main__":
    main()