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

from analysis import *
import numpy as np
import pickle
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
    ac = h5py.File('../h5ac/pils.h5', 'r')
    al = h5py.File('../h5al/pils.h5', 'r')
    frames = 10001
    Cluster = []
    freqCation, freqAnion = [], []
    for frame in range(frames):
        if frame % 100 == 0:
            print("Processing {:10.3f}".format((frame + 1) * 1.0 / frames))
        ASSOac = []
        ASSOal = []
        ASSOac = rawLoad(ac, 'htt')[frame]
        ASSOal = rawLoad(al, 'htt')[frame]
        Ncations, Nanions = [], []
        for mylithium in range(60):
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
    clusterPrint('Li_Cluster', Cluster)


if __name__ == "__main__":
    main()
