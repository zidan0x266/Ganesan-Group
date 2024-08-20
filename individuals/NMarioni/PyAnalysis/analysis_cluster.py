# This script creates a 2D array specifying all clusters based on the number of a1 and a2
# This script can be used to create custom scripts for different purposes. The 2D array is output because it can be used to gain a lot of information in post

import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
import MDAnalysis.lib.distances as distances
import numpy as np
import networkx as nx

import multiprocessing as mp
import functools
import h5py
import os
from sys import argv

script, trj_file, top_file, a1_name, a2_name, coord, t_min, t_max, step, nt = argv
coord = np.array(coord.split(), dtype = float); t_min = float(t_min); t_max = float(t_max); step = int(step); nt = int(nt)
# trj_file = .trr/.xtc file, top_file = .tpr file
# a1_name = name of reference atom to be analyzed
# a2_name = name of selection atom to be analyzed
# coord = coordination distance between a1 and a2 (Angstroms)
# t_min = start time for analysis in ps (-1 assumes the start time of the first frame)
# t_max = end time for analysis in ps (-1 assumes the end time of the last frame)
# step = frame step-size (-1 assumes a step-size of 1)
# nt = number of threads



def pairPrint(filename, pair):# print ion pair proability distribution
    with open('{}.dat'.format(filename), 'w') as anaout:
        print('# ' + filename + ' Probability', file=anaout)
        for i in range(0, len(pair)):
            print('{:10.3f} {:10.5f}'.format(i, pair[i]), file=anaout)

def printCluster2D(filename, cluster):# print ion pair proability distribution
    with open('{}.dat'.format(filename), 'w') as anaout:
        print('# Number of ion cluster with x {} and y {}'.format(a1_name, a2_name), file=anaout)
        print('# x refers to the row, starting with x = 0; y refers to the columns, starting with y = 0'.format(a1_name, a2_name), file=anaout)
        print('\n'.join(['\t'.join(['{:5.5f}'.format(cell) for cell in row]) for row in cluster]), file=anaout)



def iPair(frame):
# finds all ion clusters in a given frame
#
# Inputs: frame = time frame to be analyzed

    if frame%5000 == 0:
        print("Frame " + str(frame))

    with h5py.File('/tmp/Cluster.hdf5','r') as f:
        dset1 = f['r1']; a1 = dset1[frame]
        dset2 = f['r2']; a2 = dset2[frame]
        dset3 = f['cells']; cell = dset3[frame]

     # Clusters are defined by a1 associated with a2
    distpair = distances.capped_distance(a1, a2, coord[1], box=cell)[0]
    pair = np.zeros((len(a1), len(a2)), dtype = int)
    for index, i in enumerate(distpair[:,0]):
        pair[i,distpair[index,1]] = 1
    
     # Clusters are defined by a1 associated with a1
    #distpair = distances.self_capped_distance(a1, coord[0], box=cell)[0]
    #pair1 = np.zeros((len(a1), len(a1)), dtype = int)
    #for index, i in enumerate(distpair[:,0]):
    #    pair1[i,distpair[index,1]] = 1
    #    pair1[distpair[index,1],i] = 1

     # Clusters are defined by a2 associated with a2
    #distpair = distances.self_capped_distance(a2, coord[2], box=cell)[0]
    #pair2 = np.zeros((len(a2), len(a2)), dtype = int)
    #for index, i in enumerate(distpair[:,0]):
    #    pair2[i,distpair[index,1]] = 1
    #    pair2[distpair[index,1],i] = 1
    


    G = nx.Graph()
    for j in range(np.shape(pair)[0] + np.shape(pair)[1]):
        G.add_node(j)
    
    for j in range(np.shape(pair)[0]):
         # Uncomment if: Clusters are defined by a1 associated with a1
        #for k in range(j+1,np.shape(pair)[0]):
        #    if (pair1[j,k] == 1):
        #        G.add_edge(j, k)

        for k in range(np.shape(pair)[1]):
            if (pair[j,k] == 1):
                G.add_edge(j, np.shape(pair)[0] + k)
    
     # Uncomment if: Clusters are defined by a2 associated with a2
    #for j in range(np.shape(pair)[1]):
    #    for k in range(j+1,np.shape(pair)[1]):
    #        if (pair2[j,k] == 1):
    #            G.add_edge(np.shape(pair)[0] + j, np.shape(pair)[0] + k)


    
    n_clust = np.zeros([100,100], dtype=int)
    for z in sorted(nx.connected_components(G), key = len, reverse=True):
        z = np.array(list(z))
        i = len(z[z < np.shape(pair)[0]]); j = len(z[z >= np.shape(pair)[0]])
        n_clust[i,j] += 1
    G.clear()

    return n_clust



def load_TRR():
# loads in the trajectory and saves the necessary data to a temporary h5py file

    global t_min, t_max, step

    uta = mda.Universe(top_file, trj_file)

    a1 = uta.select_atoms("name " + a1_name)
    a2 = uta.select_atoms("name " + a2_name)

    if t_min == -1:
        t_min = uta.trajectory[0].time
    if t_max == -1:
        t_max = uta.trajectory[-1].time
    if step < 1:
        step = 1
    dt = np.round((uta.trajectory[1].time - uta.trajectory[0].time),3)
    frame_ids = np.arange(int((t_min - uta.trajectory[0].time)/dt), int((t_max - uta.trajectory[0].time)/dt + 1), step)
    dt = dt*step

    r1 = []; r2 = []; cells = []
    for frame in frame_ids:

        if frame%5000 == 0:
            print("Frame " + str(frame))
        
        ts = uta.trajectory[frame]
        cell = ts.dimensions
        cells.append(cell)

        r1.append(a1.positions)
        r2.append(a2.positions)

    with h5py.File('/tmp/Cluster.hdf5','w') as f:
        dset1 = f.create_dataset("cells", data = cells)
        dset2 = f.create_dataset("r1", data=r1)
        dset3 = f.create_dataset("r2", data=r2)
        dset4 = f.create_dataset("frames", data=frame_ids)



def main(trj_file, top_file, a1_name, a2_name, coord, t_min, t_max, step, nt):

    # Load in the trajectory file
    if not os.path.exists('/tmp/Cluster.hdf5'):
        load_TRR()
        exit()

    with h5py.File('/tmp/Cluster.hdf5','r') as f:
        dset1 = f['frames']; frame_ids = dset1[:]

    print("Pairing Analysis")
    pool = mp.Pool(processes=nt)
    func = functools.partial(iPair)
    cluster = pool.map(func, range(len(frame_ids)))
    pool.close()
    pool.join()
    cluster = np.array(cluster); cluster = np.mean(cluster, axis = 0)
    while np.all(cluster[-1,:] == 0):
        cluster = cluster[:-1,:]
    while np.all(cluster[:,-1] == 0):
        cluster = cluster[:,:-1]

     # Clust2D is an array Cxy, where x = number of a1, y = number of a2, Cxy = number of clusters with x and y
    printCluster2D('clust2D_{}_{}'.format(a1_name, a2_name), cluster)
    
    # Deletes the temporary h5py file
    os.remove('/tmp/Cluster.hdf5')

if __name__ == "__main__":
    main(trj_file, top_file, a1_name, a2_name, coord, t_min, t_max, step, nt)
