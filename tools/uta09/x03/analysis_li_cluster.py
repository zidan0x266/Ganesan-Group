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
import h5py
import networkx as nx
from collections import Counter


def clusterPrint(filename, data):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# Size', file=anaout)
        for i in range(len(data)):
            print('{:5d}'.format(data[i]), file=anaout)


def ionpairPrint(filename, xdata, ydata):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# Property Probability', file=anaout)
        for i in range(len(xdata)):
            print('{:5d} {:10.5f}'.format(xdata[i], ydata[i]), file=anaout)


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


def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def Graphs(param1, param2):
    g = nx.Graph()
    g.add_nodes_from(param1)
    g.add_edges_from(param2)
    #pos = nx.nx_pydot.graphviz_layout(g)
    g_degree = g.degree()
    #return g, pos, g_degree
    return g, g_degree


def LcenterNode(lithium, ASSOac, ASSOal):
    asso_connections = set()
    asso_activated = set()
    atom = lithium  ## lithium number lithium + 1
    canion = []
    for pair in range(len(ASSOal)):
        atom1, atom2 = ASSOal[pair]
        if atom2 == atom:
            asso_connections.add(tuple(sorted(('A' + '{}'.format(ASSOal[pair][0]), 'L' + '{}'.format(ASSOal[pair][1])))))
            asso_activated.add('A' + '{}'.format(ASSOal[pair][0]))
            asso_activated.add('L' + '{}'.format(ASSOal[pair][1]))
            canion.append(atom1)
    for i in canion:
        for pair in range(len(ASSOac)):
            atom1, atom2 = ASSOac[pair]
            if atom1 == i:
                asso_connections.add(tuple(sorted(('A' + '{}'.format(ASSOac[pair][0]), 'P' + '{}'.format(ASSOac[pair][1])))))
                asso_activated.add('A' + '{}'.format(ASSOac[pair][0]))
                asso_activated.add('P' + '{}'.format(ASSOac[pair][1]))
    return asso_activated, asso_connections


def main():
    # load ion pair data from .h5 files
    ac = h5py.File('h5ac/pils.h5', 'r')
    al = h5py.File('h5al/pils.h5', 'r')
    NLITHIUM = 400
    nrefs = 21
    Cluster = []
    freqCation, freqAnion = [], []
    nframes = 10001
    for frame in np.linspace(0, nframes - 1, nrefs, dtype = int):
        print("Currently processing to {}%!".format(frame/nframes))
        ASSOac = []
        ASSOal = []
        ASSOac = rawLoad(ac, 'htt')[frame]
        ASSOal = rawLoad(al, 'htt')[frame]
        Ncations, Nanions = [], []
        for mylithium in range(NLITHIUM):
            tcation, tanion = 0, 0
            asso_activated, asso_connections = LcenterNode(mylithium, ASSOac, ASSOal)
            #g, pos, g_degree = Graphs(asso_activated, asso_connections)
            g, g_degree = Graphs(asso_activated, asso_connections)
            Cluster.append(len(g.degree))
            for particle in g.nodes:
                if particle[0] == 'L':
                    continue
                elif particle[0] == 'P':
                    tcation += 1
                elif particle[0] == 'A':
                    tanion += 1
            Ncations.append(tcation)
            Nanions.append(tanion)
        freqCation.append(Counter(Ncations))
        freqAnion.append(Counter(Nanions))
    cationNum, cationPro = assoFreq(freqCation)
    ionpairPrint('Li_Cation', cationNum, cationPro)
    anionNum, anionPro = assoFreq(freqAnion)
    ionpairPrint('Li_Anion', anionNum, anionPro)
    clusterPrint('Licluster', Cluster)


if __name__ == "__main__":
    main()
