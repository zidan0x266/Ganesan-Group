import argparse, os
import pickle
import h5py
import networkx as nx

def bakcalcNode(h5, frame):
    cg_connections = set()
    cg_activated = set()
    for b1, b2 in cl[frame]:
        if b1 != -1 and b2 != -1:
            cg1, cg2 = map(at_cg_bead.get, (b1, b2))
            cg_connections.add(tuple(sorted((cg1, cg2))))
            cg_activated.add(cg1)
            cg_activated.add(cg2)

    return cg_activated, cg_connections


def calcNode(ASSOac):
    asso_connections = set()
    asso_activated = set()
    for i in range(len(ASSOac)):
        atom1, atom2 = ASSOac[i]
        asso_connections.add(tuple(sorted((atom1, atom2))))
        asso_activated.add(atom1)
        asso_activated.add(atom2)
    return asso_activated, asso_connections


def main():
    parser = argparse.ArgumentParser('Output data for projection FIGURE')
    parser.add_argument('--h5', help='Input original H5 file')
    parser.add_argument('--out', help='Output dat file')
    parser.add_argument('--frame', type=int, default=-1)
    args = parser.parse_args()

    h5 = h5py.File(args.h5, 'r')
    set1, set2 = calcNode(h5)
    with open(args.out, 'wb') as fp:
        pickle.dump([set1, set2], fp)

if __name__ == '__main__':
    main()