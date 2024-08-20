"""
Copyright (C) 2018-2020 Zidan Zhang <zhangzidan@gmail.com>
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

import h5py
import numpy as np
import multiprocessing as mp
import functools

def acfCtPrint(filename, aacf):
    """
    This function print the intermitent Association Lifetime
    """
    acfcout = open(filename + '_Ct.dat', 'w')
    print('# C_t', file=acfcout)
    for i in range(0, len(aacf)):
        print('{:10.5f}'.format(aacf[i]), file=acfcout)
    acfcout.close()


def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def getC_i(ASSO, t0, t):
    """
    This function give the intermitent Association Lifetime
    C_i = <h(t0)h(t)>/<h(t0)> between t0 and t
    """
    if ASSO[t0].shape[0] == 0:
        return 0.0
    
    Hit = np.intersect1d(ASSO[t0], ASSO[t])
    C_i = len(Hit) / float(len(ASSO[t0]))

    return C_i


def intervC_i(h5filename, dataset_name, t0, tf, reservoir, dt):
    """
    This function gets all the data for the h(t0)h(t0+dt), where
    t0 = 1,2,3,...,tf. This function give us a point of the final plot
    C(t) vs t,
    """
    ASSO, h5 = rawLoad(h5filename, dataset_name)
    a = 0
    count = 0
    for t0 in reservoir:
        if t0 + dt <= tf:
            a += getC_i(ASSO, t0, t0 + dt)
            count += 1
    h5.close()
    return a / count


def finalGetC_i(h5filename, dataset_name, t0, trestart, nt):
    """
    This function gets the final data of the C_i graph.
    """
    #ASP, h5 = rawLoad(h5filename1, dataset_name)
    #tf = ASP.shape[0] - 1
    tf = 5000
    #num_ASP = ASP.shape[0]
    num_ASP = 5001
    #h5.close()

    pool = mp.Pool(processes=nt)
    reservoir = [x * trestart for x in range(tf) if x * trestart < tf]
    func = functools.partial(intervC_i, h5filename, dataset_name, t0, tf, reservoir)
    return_data = pool.map(func, list(range(num_ASP)))
    return return_data


def main():
    print('calculating C(t)')
    C_t = finalGetC_i('ac.h5', 'htt', 0, 1, 48)
    acfCtPrint('subaacf', C_t)


if __name__ == "__main__":
    main()