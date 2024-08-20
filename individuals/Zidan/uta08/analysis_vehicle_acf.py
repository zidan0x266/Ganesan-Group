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

import multiprocessing as mp
import functools
import numpy as np
import h5py
import time
import MDAnalysis as mda


def acfPrint(filename, acf):
    """
    This function print the continous Association Lifetime
    """
    acfcout = open(filename + '.xvg', 'w')
    print('# Time S_t', file=acfcout)
    t = np.linspace(0, len(acf) - 1, len(acf))
    for i in range(0, len(acf)):
        print('{} {:10.5f}'.format(t[i], acf[i]), file=acfcout)
    acfcout.close()


def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def get_prob(old, new):
    if list(new) == list(old):
        prob = 1.0
    else:
        prob = 0.0
    return prob


def get_state(prev_s, temp_s):
    if prev_s == 1.0:
        state = temp_s
    else:
        state = 0.0
    return state


def getC_c(Vehicle, t0, t):
    """
    This function give the continous Association Lifetime
    C_c = <h(t0)H(t)>/<h(t0)> between t0 and t.
    """
    if t0 == t:
        return 1.0

    dt = 1
    old_sol_shell = Vehicle[t0]
    prev_state = [1.0 for x in range(Vehicle.shape[1])]

    while t0 + dt <= t:
        new_sol_shell = Vehicle[t0 + dt]
        temp_state = [get_prob(old, new) for old, new in zip(old_sol_shell, new_sol_shell)]
        curr_state = [get_state(prev_s, temp_s) for prev_s, temp_s in zip(prev_state, temp_state)]
        old_sol_shell = new_sol_shell
        prev_state = curr_state
        t0 += dt       
    C_c = np.sum(curr_state) / Vehicle.shape[1]
    if Vehicle[t0].shape[0] == 0:
        return 0
    else:
        return C_c


def intervC_c(h5filename, dataset_name, tf, reservoir, dt):
    """
    This function gets all the data for the h(t0)H(t0+dt), where
    t0 = 1,2,3,...,tf. This function give us one point of the final plot
    S(t) vs t.
    """
    Vehicle, h5 = rawLoad(h5filename, dataset_name)
    a = 0
    count = 0
    for t0 in reservoir:
        if t0 + dt <= tf:
            if t0 == t0 + dt:
                break
            a += getC_c(Vehicle, t0, t0 + dt)
            count += 1

    h5.close()

    if count == 0:
        return 1.0

    return a / count


def finalGetC_c(h5filename, dataset_name, trestart, nt):
    """
    This function gets the final data of the C_c graph.
    """
    Vehicle, h5 = rawLoad(h5filename, dataset_name)
    tf = Vehicle.shape[0] - 1
    num_Vehicle = Vehicle.shape[0]
    print(('finalGetC_c tf: {}'.format(tf)))
    h5.close()

    pool = mp.Pool(processes=nt)
    reservoir = [x * trestart for x in range(tf) if x * trestart < tf]
    func = functools.partial(intervC_c, h5filename, dataset_name, tf, reservoir)
    return_data = pool.map(func, list(range(num_Vehicle)))
    return return_data


def main():
    # generic parameters
    trestart = 100
    nt = 56  # number of processors
    h5filename = 'solvation.h5'
    dataset_name = 'shell'
    # generate the coordinates object
    S_t = finalGetC_c(h5filename, dataset_name, trestart, nt)
    acfPrint('acf_St', S_t)

if __name__ == "__main__":
    main()
