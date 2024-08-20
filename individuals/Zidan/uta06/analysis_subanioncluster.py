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
from analysis import *
import networkx as nx
from collections import Counter


def clusterPrint(filename, data):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# Size', file=anaout)
        for i in range(len(data)):
            print('{:5d}'.format(data[i]), file=anaout)


def Graphs(param1, param2):
    g = nx.Graph()
    g.add_nodes_from(param1)
    g.add_edges_from(param2)
    #pos = nx.nx_pydot.graphviz_layout(g)
    g_degree = g.degree()
    #return g, pos, g_degree
    return g, g_degree


def AcenterNode(anion, aniontype, ASSOac, ASSOal):
    asso_connections = set()
    asso_activated = set()
    atom = anion
    anions = []
    if aniontype == 1:
        for pair in ASSOac:
            atom1, atom2 = pair
            if atom1 == atom:
                asso_connections.add(tuple(sorted(('A' + '{}'.format(pair[0]), 'P' + '{}'.format(pair[1])))))
                asso_activated.add('A' + '{}'.format(pair[0]))
                asso_activated.add('P' + '{}'.format(pair[1]))
    elif aniontype == 2:
        for pair in ASSOal:
            atom1, atom2 = pair
            if atom1 == atom:
                asso_connections.add(tuple(sorted(('A' + '{}'.format(pair[0]), 'L' + '{}'.format(pair[1])))))
                asso_activated.add('A' + '{}'.format(pair[0]))
                asso_activated.add('L' + '{}'.format(pair[1]))
                anions.append(atom1)
        for i in anions:
            for pair in ASSOac:
                atom1, atom2 = pair
                if atom1 == i:
                    asso_connections.add(tuple(sorted(('A' + '{}'.format(pair[0]), 'P' + '{}'.format(pair[1])))))
                    asso_activated.add('A' + '{}'.format(pair[0]))
                    asso_activated.add('P' + '{}'.format(pair[1]))
    elif aniontype == 3:
        for pair in ASSOal:
            atom1, atom2 = pair
            if atom1 == atom:
                asso_connections.add(tuple(sorted(('A' + '{}'.format(pair[0]), 'L' + '{}'.format(pair[1])))))
                asso_activated.add('A' + '{}'.format(pair[0]))
                asso_activated.add('L' + '{}'.format(pair[1]))
    return asso_activated, asso_connections


def getgroups(Nanions, ASSOac, ASSOal):  # get the co-coordination anions
    at1, at2 = [], []
    anion = list(set([x[0] for x in ASSOal]))  # find the anions in co-coordination
    for pair in ASSOac:
        atom = pair[0]
        if atom not in anion:
            at1.append(atom)  # anions only associated with cation
        else:
            at2.append(atom)  # anions in the co-coordination
    type1 = list(set(at1))
    type2 = list(set(at2))
    anions = set(np.linspace(0, Nanions - 1, Nanions, dtype = int))
    type3 = list(anions - set(at1) - set(at2))
    return type1, type2, type3


def main():
    Nanions = 480  # number of anions
    frames = 100  # number of frames
    ac = h5py.File('../h5ac/pils.h5', 'r')
    al = h5py.File('../h5al/pils.h5', 'r')
    ASSOac = rawLoad(ac, 'htt')
    ASSOal = rawLoad(al, 'htt')
    Cluster_t1, Cluster_t2, Cluster_t3 = [], [], []
    Ncations_t2, Nlithiums_t2 = [], []
    for frame in range(frames):
        if (frame + 1) % 10 == 0:
            print("Processing {:10.3f}".format((frame + 1) * 1.0 / frames))
        type1, type2, type3 = getgroups(Nanions, ASSOac[frame], ASSOal[frame])
        for myanion in range(Nanions):
            if myanion in type1:
                tcation, tlithium = 0, 0
                asso_activated, asso_connections = AcenterNode(myanion, 1, ASSOac[frame], ASSOal[frame])
                g, g_degree = Graphs(asso_activated, asso_connections)
                Cluster_t1.append(len(g.degree))
            elif myanion in type2:
                tcation, tlithium = 0, 0
                asso_activated, asso_connections = AcenterNode(myanion, 2, ASSOac[frame], ASSOal[frame])
                g, g_degree = Graphs(asso_activated, asso_connections)
                Cluster_t2.append(len(g.degree))
                for particle in g.nodes:
                    if particle[0] == 'L':
                        tlithium += 1
                    elif particle[0] == 'P':
                        tcation += 1
                    elif particle[0] == 'A':
                        continue
                Ncations_t2.append(tcation)
                Nlithiums_t2.append(tlithium)
            elif myanion in type3:
                tcation, tlithium = 0, 0
                asso_activated, asso_connections = AcenterNode(myanion, 3, ASSOac[frame], ASSOal[frame])
                g, g_degree = Graphs(asso_activated, asso_connections)
                Cluster_t3.append(len(g.degree))
    clusterPrint('Anion_Cluster_t1', Cluster_t1)
    clusterPrint('Anion_Cluster_t2', Cluster_t2)
    clusterPrint('Anion_Cluster_t3', Cluster_t3)
    clusterPrint('Anion_t2_cations', Ncations_t2)
    clusterPrint('Anion_t2_lithiums', Nlithiums_t2)
    

if __name__ == "__main__":
    main()
