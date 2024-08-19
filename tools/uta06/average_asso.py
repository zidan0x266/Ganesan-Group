"""
Copyright (C) 2018-2021 Zidan Zhang <zhangzidan@gmail.com>

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

def outPrint(filename, xdata, ydata, yerror):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# V P(V)', file=anaout)
        for i in range(0, len(xdata)):
            print('{:10.3f} {:10.5f} {:10.5f}'.format(xdata[i], ydata[i], yerror[i]), file=anaout)


def averaging(npoints, data):
    nsamples = 10
    samples = np.zeros((nsamples, npoints))
    for counter in range(nsamples):
        for i in range(npoints):
            for j in range(len(data[counter])):
                if data[counter][:, 0][j] == i:
                    samples[counter][i] = data[counter][:, 1][j]
    return samples


def main():
    func = 0
    conc = 'c060'
    folder = ['s01', 's02', 's03', 's04', 's05',
              's06', 's07', 's08', 's09', 's10']
    if func == 0:
        maxN = []
        for dir in folder:
            data = np.loadtxt('../../' + dir + '/' + conc + '/licluster/Li_Anion_t1.xvg', usecols = 1)
            maxN.append(len(data))
        npoints = np.max(maxN) * 2
        raw = []
        for dir in folder:
            raw.append(np.loadtxt('../../' + dir + '/' + conc + '/licluster/Li_Anion_t1.xvg'))
        samples = averaging(npoints, raw)
        xdata = np.linspace(0, npoints - 1, npoints)
        ydata = np.average(samples, axis = 0)
        yerror = np.std(samples, axis = 0)
        outPrint('Li_Anion_t1_' + conc, xdata, ydata, yerror)
        func = 1
    if func == 1:
        maxN = []
        for dir in folder:
            data = np.loadtxt('../../' + dir + '/' + conc + '/licluster/Li_Anion_t2.xvg', usecols = 1)
            maxN.append(len(data))
        npoints = np.max(maxN) * 2
        raw = []
        for dir in folder:
            raw.append(np.loadtxt('../../' + dir + '/' + conc + '/licluster/Li_Anion_t2.xvg'))
        samples = averaging(npoints, raw)
        xdata = np.linspace(0, npoints - 1, npoints)
        ydata = np.average(samples, axis = 0)
        yerror = np.std(samples, axis = 0)
        outPrint('Li_Anion_t2_' + conc, xdata, ydata, yerror)
        func = 2
    if func == 2:
        maxN = []
        for dir in folder:
            data = np.loadtxt('../../' + dir + '/' + conc + '/licluster/Li_Cation_t1.xvg', usecols = 1)
            maxN.append(len(data))
        npoints = np.max(maxN)
        raw = []
        for dir in folder:
            raw.append(np.loadtxt('../../' + dir + '/' + conc + '/licluster/Li_Cation_t1.xvg'))
        samples = averaging(npoints, raw)
        xdata = np.linspace(0, npoints - 1, npoints)
        ydata = np.average(samples, axis = 0)
        yerror = np.std(samples, axis = 0)
        outPrint('Li_Cation_t1_' + conc, xdata, ydata, yerror)
        func = 3


if __name__ == "__main__":
    main()