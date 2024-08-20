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

def outPrint(filename, xdata, ydata, yerror):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# V P(V)', file=anaout)
        for i in range(0, len(xdata)):
            print('{:10.3f} {:10.5f} {:10.5f}'.format(xdata[i], ydata[i], yerror[i]), file=anaout)


def gathering(npoints, data):
    nsamples = 5
    samples = np.zeros((nsamples, npoints))
    for counter in range(nsamples):
        for i in range(npoints):
            for j in range(len(data[counter])):
                if data[counter][:, 0][j] == i:
                    samples[counter][i] = data[counter][:, 1][j]
    return samples


def averaging(data, npoints):
    ydata, yerror = [], []
    for i in range(npoints):
        tmp = [x for x in data[:,i] if x != 0]
        if len(tmp) == 0:
            ydata.append(0.0)
            yerror.append(0.0)
        else:
            ydata.append(np.average(tmp))
            yerror.append(np.std(tmp))
    return ydata, yerror


def main():
    func = 0  # 0 for chain association and 1 for pair assiciation
    xmax = 25
    joblist = ['li_anion', 'li_cation', 'li_chain', 'li_coco']
    folder = ['s01', 's02', 's03', 's04', 's05']
    for job in joblist:
        raw = []
        for dir in folder:
            raw.append(np.loadtxt('../' + dir + '/CONCENTRATION/linformation/' + job + '.dat'))
        samples = gathering(xmax, raw)
        xdata = np.linspace(0, xmax - 1, xmax)
        ydata, yerror = averaging(samples, xmax)
        outPrint(job + '_CONCENTRATION', xdata, ydata, yerror)


if __name__ == "__main__":
    main()