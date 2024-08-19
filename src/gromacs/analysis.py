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

import multiprocessing as mp

import functools
import h5py
import numpy as np


def assoFreq(frequence):
    """
    This function gets the normolized probability for association events.
    """
    fNum = frequence[0]
    for i in range(1, len(frequence)):
        fNum += frequence[i]
    numAsso = {}
    for k, v in fNum.items():
        numAsso[k] = v / float(len(frequence))
    totalAsso = sum(numAsso.values())
    avgAsso = {}
    for k, v in numAsso.items():
        avgAsso[k] = v / float(totalAsso)
    pairNum = []
    pairPro = []
    for k in sorted(avgAsso.keys()):
        pairNum.append(k)
        pairPro.append(avgAsso[k])
    return pairNum, pairPro


def ds_append(ds, val):
    ds.resize(ds.shape[0]+1, axis=0)
    ds[-1] = val


def get_or_create_dataset(h5, dataset_name, dtype):
    if dataset_name not in h5:
        return h5.create_dataset(dataset_name, (0, ), maxshape=(None, ), dtype=dtype, chunks=(100, ))
    return h5[dataset_name]


def rawSave(h5, dataset_name, raw_data):
    dtype = h5py.special_dtype(vlen=np.dtype('int, int'))
    if dataset_name not in h5:
        dataset = h5.create_dataset(dataset_name, (len(raw_data), ), maxshape=(None, ), dtype=dtype, chunks=(len(raw_data), ))
    else:
        dataset = h5[dataset_name]

    for i, v in enumerate(raw_data):
        dataset[i] = v


def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def rdfPrint(filename, nbins, rdf):  # print rdf
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# r g(r)', file=anaout)
        for i in range(0, len(nbins)):
            print('{:10.3f} {:10.5f}'.format(nbins[i] / 10, rdf[i]), file=anaout)


def avgPrint(filename, xdata, ydata, yerror):  # print averaged data
    if xdata[0] != 0:
        xdata.insert(0, 0)
        ydata.insert(0, 0)
        yerror.insert(0, 0)
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# number average std', file=anaout)
        for i in range(0, len(xdata)):
            print('{:d} {:10.5e} {:10.5e}'.format(int(xdata[i]), ydata[i], yerror[i]), file=anaout)


def anaPrint(filename, arrayNum, arrayPro):  # print association properties
    anaout = open(filename + '.dat', 'w')
    print('# ' + filename + ' Probability', file=anaout)
    for i in range(0, len(arrayNum)):
        print('{} {:10.5f}'.format(arrayNum[i], arrayPro[i]), file=anaout)
    anaout.close()


def interFPrint(filename, average, std):  # print number of ions near the interface
    anaout = open(filename + '.dat', 'w')
    print('# ' + filename + ' Probability', file=anaout)
    print('{} {:10.5f}'.format(average, std), file=anaout)
    anaout.close()


def hopesPrint(filename, intramolecular, intermolecular):  # print association properties
    hopout = open(filename + '.dat', 'w')
    print('# ' + 'intramolecular' + ' intermolecular', file=hopout)
    for i in range(0, len(intramolecular)):
        print('{:10.5f} {:10.5f}'.format(intramolecular[i], intermolecular[i]), file=hopout)
    hopout.close()


def acfCtPrint(filename, aacf):
    """
    This function print the intermitent Association Lifetime
    """
    acfcout = open(filename + '_Ct.dat', 'w')
    print('# Time C_t', file=acfcout)
    tc = np.linspace(0, len(aacf) - 1, len(aacf))
    for i in range(0, len(aacf)):
        print('{} {:10.5f}'.format(tc[i] * 2, aacf[i]), file=acfcout)
    acfcout.close()


def acfStPrint(filename, aacf):
    """
    This function print the continous Association Lifetime
    """
    acfcout = open(filename + '_St.dat', 'w')
    print('# Time S_t', file=acfcout)
    tc = np.linspace(0, len(aacf) - 1, len(aacf))
    for i in range(0, len(aacf)):
        print('{} {:10.5f}'.format(tc[i] * 2, aacf[i]), file=acfcout)
    acfcout.close()


def getC_i(ASP, t0, t):
    """
    This function give the intermitent Association Lifetime
    C_i = <h(t0)h(t)>/<h(t0)> between t0 and t
    """
    if ASP[t0].shape[0] == 0:
        return 0.0

    Hit = np.intersect1d(ASP[t0], ASP[t])
    C_i = len(Hit) / float(len(ASP[t0]))

    return C_i


def intervC_i(h5filename, dataset_name, t0, tf, reservoir, dt):
    """
    This function gets all the data for the h(t0)h(t0+dt), where
    t0 = 1,2,3,...,tf. This function give us a point of the final plot
    C(t) vs t,
    """
    ASP, h5 = rawLoad(h5filename, dataset_name)
    a = 0
    count = 0
    for t0 in reservoir:
        if t0 + dt <= tf:
            a += getC_i(ASP, t0, t0 + dt)
            count += 1
    h5.close()
    return a / count


def finalGetC_i(h5filename, dataset_name, t0, trestart, nt):
    """
    This function gets the final data of the C_i graph.
    """
    ASP, h5 = rawLoad(h5filename, dataset_name)
    tf = ASP.shape[0] - 1
    num_ASP = ASP.shape[0]
    h5.close()

    pool = mp.Pool(processes=nt)
    reservoir = [x * trestart for x in range(tf) if x * trestart < tf]
    func = functools.partial(intervC_i, h5filename, dataset_name, t0, tf, reservoir)
    return_data = pool.map(func, list(range(num_ASP)))
    return return_data


def getC_c(ASP, t0, t):
    """
    This function give the continous Association Lifetime
    C_c = <h(t0)H(t)>/<h(t0)> between t0 and t.
    """
    if t0 == t:
        return 1.0

    dt = 1
    ht0 = ASP[t0]
    htt_old = ASP[t0]

    while t0 + dt <= t:
        htt = ASP[t0 + dt]
        Htt = np.intersect1d(htt_old, htt)
        # Htt = set(htt_old) & set(htt)  # comparison between current frame and previous frame
        t0 += dt
        htt_old = Htt
    C_c = len(Htt) / float(len(ht0))
    if ASP[t0].shape[0] == 0:
        return 0
    else:
        return C_c


def intervC_c(h5filename, dataset_name, t0, tf, reservoir, dt):
    """
    This function gets all the data for the h(t0)H(t0+dt), where
    t0 = 1,2,3,...,tf. This function give us one point of the final plot
    S(t) vs t.
    """
    ASP, h5 = rawLoad(h5filename, dataset_name)
    a = 0
    count = 0
    for t0 in reservoir:
        if t0 + dt <= tf:
            if t0 == t0 + dt:
                break
            a += getC_c(ASP, t0, t0 + dt)
            count += 1

    h5.close()

    if count == 0:
        return 1.0

    return a / count


def finalGetC_c(h5filename, dataset_name, t0, trestart, nt):
    """
    This function gets the final data of the C_c graph.
    """
    ASP, h5 = rawLoad(h5filename, dataset_name)
    tf = ASP.shape[0] - 1
    num_ASP = ASP.shape[0]
    print(('finalGetC_c tf: {}'.format(tf)))
    h5.close()

    pool = mp.Pool(processes=nt)
    reservoir = [x * trestart for x in range(tf) if x * trestart < tf]
    func = functools.partial(intervC_c, h5filename, dataset_name, t0, tf, reservoir)
    return_data = pool.map(func, list(range(num_ASP)))
    return return_data
