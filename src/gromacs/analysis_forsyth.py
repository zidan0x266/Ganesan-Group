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

def construction(h5, dataset_name, length, dtype):
    dtype = h5py.vlen_dtype(np.dtype('int32'))
    if dataset_name not in h5:
        #return h5.create_dataset(dataset_name, (0, length), dtype=dtype)
        return h5.create_dataset(dataset_name, (0, length), maxshape=(None, length), dtype=dtype, chunks=(100, length))
    return h5[dataset_name]


def ds_append(ds, val):
    ds.resize(ds.shape[0]+1, axis=0)
    ds[-1] = val


def rawSave(h5, dataset_name, raw_data, length):
    dtype = h5py.vlen_dtype(np.dtype('int32'))
    if dataset_name not in h5:
        dataset = h5.create_dataset(dataset_name, (len(raw_data), length), dtype=dtype)
    else:
        dataset = h5[dataset_name]
    for i, v in enumerate(raw_data):
        dataset[i] = v


def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def acfCtPrint(filename, aacf):
    """
    This function print the intermitent Association Lifetime
    """
    acfcout = open(filename + '_Ct.dat', 'w')
    print('# Time X_t', file=acfcout)
    tc = np.linspace(0, len(aacf) - 1, len(aacf))
    for i in range(0, len(aacf)):
        print('{} {:10.5f}'.format(tc[i], aacf[i]), file=acfcout)
    acfcout.close()


def getC_i1(Vehicle, t0, t):
    """
    This function give the intermitent Association Lifetime
    C_i = <h(t0)h(t)>/<h(t0)> between t0 and t
    """
    if Vehicle[t0].shape[0] == 0:
        return 0.0
    
    probability = []
    for plate in range(len(Vehicle[t0])):
        overlap = np.intersect1d(Vehicle[t0][plate], Vehicle[t][plate])
        if len(Vehicle[t0][plate]) < 3:
            continue
        elif len(overlap) == len(Vehicle[t0][plate]) and len(Vehicle[t0][plate]) >= 3:
            prob = 1.0
        elif len(overlap) <= 1 and len(Vehicle[t0][plate]) >= 3:
            prob = 0.0
        else:
            prob = 0.9
        probability.append(prob)
    
    C_i = np.average(probability)

    return C_i


def get_prob(a, b):
    if b < 3:
        prob = None
    elif a == b and b >= 3:
        prob = 1.0
    elif a <= 1 and b >= 3:
        prob = 0.0
    else:
        prob = 0.9
    return prob


def getC_i2(Vehicle, t0, t):
    """
    This function give the intermitent Association Lifetime
    C_i = <h(t0)h(t)>/<h(t0)> between t0 and t
    """
    if Vehicle[t0].shape[0] == 0:
        return 0.0

    origin = [len(x) for x in Vehicle[t0]]
    overlap = [len(np.intersect1d(x, y)) for x, y in zip(Vehicle[t0], Vehicle[t])]
    temp = [get_prob(x, y) for x, y in zip(overlap, origin)]
    probability = [x for x in temp if x != None]
    C_i = np.average(probability)

    return C_i


def intervC_i(h5filename, dataset_name, t0, tf, reservoir, dt):
    """
    This function gets all the data for the h(t0)h(t0+dt), where
    t0 = 1,2,3,...,tf. This function give us a point of the final plot
    C(t) vs t,
    """
    Vehicle, h5 = rawLoad(h5filename, dataset_name)
    a = 0
    count = 0
    for t0 in reservoir:
        if t0 + dt <= tf:
            a += getC_i2(Vehicle, t0, t0 + dt)
            count += 1
    h5.close()
    return a / count


def finalGetC_i(h5filename, dataset_name, t0, trestart, nt):
    """
    This function gets the final data of the C_i graph.
    """
    Vehicle, h5 = rawLoad(h5filename, dataset_name)
    tf = Vehicle.shape[0] - 1
    num_Vehicle = Vehicle.shape[0]
    h5.close()

    pool = mp.Pool(processes=nt)
    reservoir = [x * trestart for x in range(tf) if x * trestart < tf]
    func = functools.partial(intervC_i, h5filename, dataset_name, t0, tf, reservoir)
    return_data = pool.map(func, list(range(num_Vehicle)))
    return return_data


def get_vehicle(lithium, ASSOac, ASSOal):
    asso_connections = set()
    asso_activated = set()
    atom = lithium
    inner_vehicle = []  # first solvation shell
    intra_vehicle = []  # second solvation shell
    linker = []
    for pair in range(len(ASSOal)):
        atom1, atom2 = ASSOal[pair]
        if atom2 == atom:
            inner_vehicle.append(atom1)
            linker.append(atom1)
    for pair in range(len(ASSOac)):
        atom1, atom2 = ASSOac[pair]
        if atom1 in linker:
            intra_vehicle.append(atom2)
    return list(set(inner_vehicle)), list(set(intra_vehicle))


def vehicular(Nlithiums, frame_id):
    H5ac = h5py.File('../h5ac/pils.h5', 'r')
    H5al = h5py.File('../h5al/pils.h5', 'r')
    ASSOac = rawLoad(H5ac, 'htt')
    ASSOal = rawLoad(H5al, 'htt')
    Inner_vehicle, Intra_vehicle = [], []
    for mylithium in range(Nlithiums):
        inner_vehicle, intra_vehicle = get_vehicle(mylithium, ASSOac[frame_id], ASSOal[frame_id])
        Inner_vehicle.append(inner_vehicle)
        Intra_vehicle.append(intra_vehicle)
    return Inner_vehicle, Intra_vehicle


def main():
    Nlithiums = 540  # number of anions
    frames = 10001  # number of frames
    nt = 56  # number of processors
    trestart = 10  # frequency of the reference point
    dtype = h5py.vlen_dtype(np.dtype('int32'))
    h5 = h5py.File('vehicle.h5', 'w')  # save vehicle information
    Inner_dataset = construction(h5, 'Inner_v', Nlithiums, dtype)
    Intra_dataset = construction(h5, 'Intra_v', Nlithiums, dtype)
    #ac = h5py.File('h5ac/pils.h5', 'r')
    #al = h5py.File('h5al/pils.h5', 'r')
    #ASSOac = rawLoad(ac, 'htt')
    #ASSOal = rawLoad(al, 'htt')
    ## TODO
    #frame_ids = [frame for frame in range(frames) if (not begin or frame >= begin) and (not end or frame <= end)]
    frame_ids = [frame for frame in range(frames)]
    vehicle_analyse = functools.partial(vehicular, Nlithiums)
    pool = mp.Pool(processes=nt)
    data = pool.imap(vehicle_analyse, frame_ids)
    for Inner_v, Intra_v in data:
        ds_append(Inner_dataset, Inner_v)
        ds_append(Intra_dataset, Intra_v)
    #Inner_v, Intra_v = [], []
    #for frame in range(frames):
    #    if (frame + 1) % 10 == 0:
    #        print("Processing {:10.3f}".format(frame * 1.0 / frames))
    #    Inner_vehicle, Intra_vehicle = vehicular(Nlithiums, ASSOac, ASSOal, frame)
    #    Inner_v.append(Inner_vehicle)
    #    Intra_v.append(Intra_vehicle)
    #rawSave(h5, 'Inner_v', Inner_v, Nlithiums)
    #rawSave(h5, 'Intra_v', Intra_v, Nlithiums)
    h5.close()
    C_t = finalGetC_i('vehicle.h5', 'Inner_v', 0, trestart, nt)
    acfCtPrint('Inner_v', C_t)
    C_t = finalGetC_i('vehicle.h5', 'Intra_v', 0, trestart, nt)
    acfCtPrint('Intra_v', C_t)
    

if __name__ == "__main__":
    main()
