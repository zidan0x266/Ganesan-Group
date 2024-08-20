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


def getSlope(time, msd, msdrange):
    #params = stats.linregress(MSD[msdrange[0]:msdrange[1], 0], MSD[msdrange[0]:msdrange[1], axis])  # axis: 1 - anion; 2 - cation
    #print("Diffusivity is {:8.5e}".format(params[0] * 1e-2 / 6))
    params = stats.linregress(np.log(time[msdrange[0]:msdrange[1]]), np.log(msd[msdrange[0]:msdrange[1]]))
    return params[0]


def polynom(xinput, yinput, fitrange):
    xdata = xinput
    ydata = yinput
    time = [1 / x for x in xdata[1:]]
    MSDdt = [x / (2 * y) for x, y in zip(ydata[1:], xdata[1:])]
    msdrange = fitrange
    x = time[msdrange[0]:msdrange[1]]
    y = MSDdt[msdrange[0]:msdrange[1]]
    popt0 = np.polyfit(x, y, 0)
    popt1 = np.polyfit(x, y, 1)
    slope = getSlope(xdata[1:], ydata[1:], msdrange)
    print("{:8.5e}, {:8.5e}, {:8.5e}, {:8.5e}".format(popt0[-1] * 1e-2, popt1[-1] * 1e-2, popt1[-2], slope))
    return [popt0[-1] * 1e-2, popt1[-1] * 1e-2, popt1[-2], slope]


def polynom3D(xinput, yinput, fitrange):
    xdata = xinput
    ydata = yinput
    time = [1 / x for x in xdata[1:]]
    MSDdt = [x / (6 * y) for x, y in zip(ydata[1:], xdata[1:])]
    msdrange = fitrange
    x = time[msdrange[0]:msdrange[1]]
    y = MSDdt[msdrange[0]:msdrange[1]]
    popt0 = np.polyfit(x, y, 0)
    popt1 = np.polyfit(x, y, 1)
    slope = getSlope(xdata, ydata, msdrange)
    print("{:8.5e}, {:8.5e}, {:8.5e}, {:8.5e}".format(popt0[-1] * 1e-2, popt1[-1] * 1e-2, popt1[-2], slope))
    return [popt0[-1] * 1e-2, popt1[-1] * 1e-2, popt1[-2], slope]


def getDtensor(MSD, fitrange):
    tensorD = np.zeros((3, 3))
    #dxx = np.average(MSD[fitrange[0]:fitrange[1], 2] / (MSD[fitrange[0]:fitrange[1], 0] / 1000)) * 1e-5 / 2
    #dyy = np.average(MSD[fitrange[0]:fitrange[1], 3] / (MSD[fitrange[0]:fitrange[1], 0] / 1000)) * 1e-5 / 2
    #dzz = np.average(MSD[fitrange[0]:fitrange[1], 4] / (MSD[fitrange[0]:fitrange[1], 0] / 1000)) * 1e-5 / 2
    #dyx = np.average(MSD[fitrange[0]:fitrange[1], 5] / (MSD[fitrange[0]:fitrange[1], 0] / 1000)) * 1e-5 / 2
    #dzx = np.average(MSD[fitrange[0]:fitrange[1], 6] / (MSD[fitrange[0]:fitrange[1], 0] / 1000)) * 1e-5 / 2
    #dzy = np.average(MSD[fitrange[0]:fitrange[1], 7] / (MSD[fitrange[0]:fitrange[1], 0] / 1000)) * 1e-5 / 2
    dxx = polynom(MSD[:, 0], MSD[:, 2], fitrange)[1]
    dyy = polynom(MSD[:, 0], MSD[:, 3], fitrange)[1]
    dzz = polynom(MSD[:, 0], MSD[:, 4], fitrange)[1]
    dyx = polynom(MSD[:, 0], MSD[:, 5], fitrange)[1]
    dzx = polynom(MSD[:, 0], MSD[:, 6], fitrange)[1]
    dzy = polynom(MSD[:, 0], MSD[:, 7], fitrange)[1]
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
    file_prefix = ['msd_a', 'msd_c']
    func_type = {0: 'HO', 1: 'LA'}
    func_key = 1
    fitrange = [30000, 280000]    
    if func_type[func_key] == 'HO':
        polynoDa = []
        polynoDc = []
        rawTF = loadrawMSD(filepath, file_prefix[0])
        rawIM = loadrawMSD(filepath, file_prefix[1])
        for i in range(len(rawTF)):
            polynoDa.append(polynom3D(rawTF[i][:, 0], rawTF[i][:, 1], fitrange)[1])
            polynoDc.append(polynom3D(rawIM[i][:, 0], rawIM[i][:, 1], fitrange)[1])
        print("{:8.5e}".format(np.average(polynoDa)))
        print("{:8.5e}".format(np.average(polynoDc)))
    elif func_type[func_key] == 'LA':
        polynoDa = []
        polynoDc = []
        rawTF = loadrawMSD(filepath, file_prefix[0])
        rawIM = loadrawMSD(filepath, file_prefix[1])
        for i in range(len(rawTF)):
            tmp_vals_a, _ = la.eig(getDtensor(rawTF[i], fitrange))
            polynoDa.append(np.average([x for x in tmp_vals_a.real if x != np.min(tmp_vals_a.real)]))
            tmp_vals_c, _ = la.eig(getDtensor(rawIM[i], fitrange))
            polynoDc.append(np.average([x for x in tmp_vals_c.real if x != np.min(tmp_vals_c.real)]))
        print("{:8.5e}".format(np.average(polynoDa)))
        print("{:8.5e}".format(np.average(polynoDc)))


if __name__ == "__main__":
    main()