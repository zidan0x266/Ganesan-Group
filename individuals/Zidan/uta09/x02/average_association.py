"""
Copyright (C) 2018-2020 Zidan Zhang <zhangzidan@gmail.com>

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
    nsamples = 5
    samples = np.zeros((nsamples, npoints))
    for counter in range(nsamples):
        for i in range(npoints):
            for j in range(len(data[counter])):
                if data[counter][:, 0][j] == i:
                    samples[counter][i] = data[counter][:, 1][j]
    return samples


def processing(folder_list, data_type):
    maxN, raw = [], []
    for dir in folder_list:
        data = np.loadtxt('../../' + dir + '/domains/anion_cation_' + data_type + '.xvg', usecols = 1)
        maxN.append(len(data))
    npoints = np.max(maxN)
    for dir in folder_list:
        raw.append(np.loadtxt('../../' + dir + '/domains/anion_cation_' + data_type + '.xvg'))
    samples = averaging(npoints, raw)
    xdata = np.linspace(0, npoints - 1, npoints)
    ydata = np.average(samples, axis = 0)
    yerror = np.std(samples, axis = 0)
    outPrint('asso_' + data_type, xdata, ydata, yerror)

def main():
    funcs = [True, True]
    folder = ['01', '02', '03', '04', '05']
    chain_list = ['Chain_to', 'Chain_bk', 'Chain_if']
    pair_list = ['Pair_to', 'Pair_bk', 'Pair_if']
    if funcs[0]:
        for job in chain_list:
            processing(folder, job)
    if funcs[1]:
        for job in pair_list:
            processing(folder, job)


if __name__ == "__main__":
    main()