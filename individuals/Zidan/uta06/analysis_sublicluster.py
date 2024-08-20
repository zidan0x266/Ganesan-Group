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
import networkx as nx
from collections import Counter


def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def groupPrint(filename, anion1, anion2):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# polyIL Both', file=anaout)
        for i in range(len(anion1)):
            print('{:5d} {:5d}'.format(anion1[i], anion2[i]), file=anaout)


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


def ligroups(frame, ASSOac, ASSOal):
    at1, at2 = [], []
    coion = list(set([x[0] for x in ASSOac[frame]]))
    firstli = True
    for pair in ASSOal[frame]:
        lithium, anion = pair
        if firstli:
            liasso = []
            cur_li = lithium
            liasso.append(anion)
            firstli = False
        else:
            if lithium == cur_li:
                liasso.append(anion)
            else:
                if set(liasso) & set(coion):
                    at1.append(cur_li)
                else:
                    at2.append(cur_li)
                cur_li = lithium
                liasso = []
                liasso.append(anion)
        if pair == ASSOal[frame][-1]:
            if set(liasso) & set(coion):
                at1.append(cur_li)
            else:
                at2.append(cur_li)
    return list(set(at1)), list(set(at2))


def main():
    ac = h5py.File('../h5ac/pils.h5', 'r')
    al = h5py.File('../h5al/pils.h5', 'r')
    la = h5py.File('../h5la/pils.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')
    ASSOla = rawLoad(la, 'htt')
    polyIL, Only = [], []
    frames = 10001
    freqCation_t1, freqAnion_t1 = [], []
    freqCation_t2, freqAnion_t2 = [], []
    Cluster_t1, Cluster_t2, Cluster_t3 = [], [], []
    for time in range(frames):
        if (time + 1) % 100 == 0:
            print("Processing {:10.3f}".format((time) * 1.0 / frames))
        poly, only = ligroups(time, ASSOac, ASSOla)
        polyIL.append(len(poly))
        Only.append(len(only))
        Ncations_t1, Nanions_t1 = [], []
        Ncations_t2, Nanions_t2 = [], []
        for mylithium in range(180):
            if mylithium in poly:
                tcation, tanion = 0, 0
                asso_activated, asso_connections = LcenterNode(mylithium, ASSOac[time], ASSOal[time])
                g, g_degree = Graphs(asso_activated, asso_connections)
                Cluster_t1.append(len(g.degree))
                for particle in g.nodes:
                    if particle[0] == 'L':
                        continue
                    elif particle[0] == 'P':
                        tcation += 1
                    elif particle[0] == 'A':
                        tanion += 1
                Ncations_t1.append(tcation)
                Cluster_t3.append(tanion + 1)  # with the mylithium
                Nanions_t1.append(tanion)
            elif mylithium in only:
                tcation, tanion = 0, 0
                asso_activated, asso_connections = LcenterNode(mylithium, ASSOac[time], ASSOal[time])
                g, g_degree = Graphs(asso_activated, asso_connections)
                Cluster_t2.append(len(g.degree))
                for particle in g.nodes:
                    if particle[0] == 'L':
                        continue
                    elif particle[0] == 'P':
                        tcation += 1
                    elif particle[0] == 'A':
                        tanion += 1
                Ncations_t2.append(tcation)
                Nanions_t2.append(tanion)
        if len(Ncations_t1) != 0:
            freqCation_t1.append(Counter(Ncations_t1))
        if len(Nanions_t1) != 0:
            freqAnion_t1.append(Counter(Nanions_t1))
        if len(Ncations_t2) != 0:
            freqCation_t2.append(Counter(Ncations_t2))
        if len(Nanions_t2) != 0:
            freqAnion_t2.append(Counter(Nanions_t2))
    cationNum_t1, cationPro_t1 = assoFreq(freqCation_t1)
    ionpairPrint('Li_Cation_t1', cationNum_t1, cationPro_t1)
    cationNum_t2, cationPro_t2 = assoFreq(freqCation_t2)
    ionpairPrint('Li_Cation_t2', cationNum_t2, cationPro_t2)
    anionNum_t1, anionPro_t1 = assoFreq(freqAnion_t1)
    ionpairPrint('Li_Anion_t1', anionNum_t1, anionPro_t1)
    anionNum_t2, anionPro_t2 = assoFreq(freqAnion_t2)
    ionpairPrint('Li_Anion_t2', anionNum_t2, anionPro_t2)
    clusterPrint('Li_Cluster_t1', Cluster_t1)
    clusterPrint('Li_Cluster_t2', Cluster_t2)
    clusterPrint('Li_Cluster_t3', Cluster_t3)


if __name__ == "__main__":
    main()
