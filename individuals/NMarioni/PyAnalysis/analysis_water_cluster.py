# This scripts performs a water cluster analysis determined the cluster size and probability distributions, as well as the fraction of water molecules that are percolated in at least 1 dimension

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

script, trj_file, top_file, aw_name, coord, t_min, t_max, step, nt = argv
coord = float(coord); t_min = float(t_min); t_max = float(t_max); step = int(step); nt = int(nt)
# trj_file = .trr/.xtc file, top_file = .tpr file
# aw_name = name of reference atom to be analyzed
# NOTE: If you want the MSD of the CoM of a molecule, send in aX_name as 'resname MOL',
#       where MOL is the name of the residue of the molecule of interest
# coord = coordination distance between waters (Angstroms)
# t_min = start time for analysis in ps (-1 assumes the start time of the first frame)
# t_max = end time for analysis in ps (-1 assumes the end time of the last frame)
# step = frame step-size (-1 assumes a step-size of 1)
# nt = number of threads



def pairPrint(filename, x, arr):# print ion pair proability distribution
    with open('{}.dat'.format(filename), 'w') as anaout:
        print('# Clust_Size Number', file=anaout)
        for i in range(0, len(arr)):
            print('{:10.3f} {:10.5f}'.format(x[i], arr[i]), file=anaout)

def pairPrint2(filename, arr):# print ion pair proability distribution
    with open('{}.dat'.format(filename), 'w') as anaout:
        print('# Clust_Size Number', file=anaout)
        for i in range(0, len(arr)):
            print('{:10.3f} {:10.5f}'.format(i, arr[i]), file=anaout)

def pairPrint3(filename, arr):# print ion pair proability distribution
    with open('{}.dat'.format(filename), 'w') as anaout:
        print('# N R(N)', file=anaout)
        for i in range(0, len(arr)):
            print('{:10.3f} {:10.5f}'.format(i, arr[i]), file=anaout)



def iPair(frame):
# finds all ion clusters in a given frame
#
# Inputs: frame = time frame to be analyzed

    if frame%5000 == 0:
        print("Frame " + str(frame))

    with h5py.File('/tmp/Cluster.hdf5','r') as f:
        dset1 = f['r1']; a1 = dset1[frame]
        dset3 = f['cells']; cell = dset3[frame]
    
     # Clusters are defined by a1 associated with a1
    distpair = distances.self_capped_distance(a1, coord, box=cell)[0]
    pair = np.zeros((len(a1), len(a1)), dtype = int)
    for index, i in enumerate(distpair[:,0]):
        pair[i,distpair[index,1]] = 1
        pair[distpair[index,1],i] = 1
    


    G = nx.Graph()
    G.add_nodes_from(range(len(a1)))

    x,y = np.where(pair==1)
    G.add_edges_from(np.array([x,y]).T)
    #G = nx.Graph()
    #for j in range(len(a1)):
    #    G.add_node(j)
    #
    #for j in range(len(a1)):
    #     # Uncomment if: Clusters are defined by a1 associated with a1
    #    for k in range(j+1, len(a1)):
    #        if (pair[j,k] == 1):
    #            G.add_edge(j, k)



    clusters = np.zeros([10000])
    max_dist = np.zeros([10000])
    count = np.zeros([10000])

    max_clust = len(sorted(nx.connected_components(G), key = len, reverse=True)[0])

    P_num = 0; P_denom = 0
    for z in sorted(nx.connected_components(G), key = len, reverse=True):
        z = np.array(list(z))
        clusters[len(z)] += 1

        if len(z) > 1:
            max_dist[len(z)] += np.amax(distances.self_distance_array(a1[z], box=cell))/10/2
            count[len(z)] += 1

        if len(z) >= np.floor(cell[0]/coord):
            Percolation = Check_Percolation(a1[z], cell)
            if Percolation == 1:
                P_num += len(z)
        P_denom += len(z)
    P_frac = (np.zeros([10000]) + 1) * (P_num / P_denom)
        
    clust = clusters / np.sum(clusters); S_num = 0.0; S_denom = 0.0
    for z in sorted(nx.connected_components(G), key = len, reverse=True):
        if len(z) < max_clust:
            S_num += (clust[len(z)]*len(z)*len(z))
            S_denom += (clust[len(z)]*len(z))
    if S_denom == 0:
        S_max = (np.zeros([10000]) + 1) * 0
    else:
        S_max = (np.zeros([10000]) + 1) * (S_num / S_denom)
        
    G.clear()

    return [clusters, max_dist, count, S_max, P_frac]



def Check_Percolation(arr, cell):

    arr_0 = np.array(arr)
     # Clusters are defined by a1 associated with a1
    distpair = distances.self_capped_distance(arr_0, coord)[0]
    pair = np.zeros((len(arr_0), len(arr_0)), dtype = int)
    for index, j in enumerate(distpair[:,0]):
        pair[j,distpair[index,1]] = 1
        pair[distpair[index,1],j] = 1

    H = nx.Graph()
    for j in range(len(arr_0)):
        H.add_node(j)
    for j in range(len(arr_0)):
         # Uncomment if: Clusters are defined by a1 associated with a1
        for k in range(j+1, len(arr_0)):
            if (pair[j,k] == 1):
                H.add_edge(j, k)

    for z in sorted(nx.connected_components(H), key = len, reverse=True):
        z = np.array(list(z))
        clust = z
        break


    dim = np.zeros((3,2,100)) - 1
    for i in range(3):
        arr_1 = np.array(arr_0); arr_1[:,i] += cell[i]
        arr_2 = np.array(arr_0); arr_2[:,i] += (2*cell[i])
        
        distpair0 = distances.capped_distance(arr_1[clust], arr_0, coord)[0][:,1]
        distpair2 = distances.capped_distance(arr_1[clust], arr_2, coord)[0][:,1]

        c0 = 0; c2 = 0
        for k, z in enumerate(sorted(nx.connected_components(H), key = len, reverse=True)):
            z = np.array(list(z))
            for j in distpair0:
                if j in z and not (k in dim[i,0,:]):
                    dim[i,0,c0] = k
                    c0 += 1

            for j in distpair2:
                if j in z and not (k in dim[i,1,:]):
                    dim[i,1,c2] = k
                    c2 += 1


    for i in range(6):
        d1 = dim[int(i/2.0),i%2,:]; d1 = d1[d1 != -1]
        if len(d1) == 0:
            continue
        for j in range(i+1,6):
            d2 = dim[int(j/2.0),j%2,:]; d2 = d2[d2 != -1]
            if len(d2) == 0:
                continue

            for k in d1:
                if k in d2:
                    return 1
    
    return 0



def load_TRR():
# loads in the trajectory and saves the necessary data to a temporary h5py file

    global t_min, t_max, step

    uta = mda.Universe(top_file, trj_file)

    a1 = uta.select_atoms("name " + aw_name)

    if t_min == -1:
        t_min = uta.trajectory[0].time
    if t_max == -1:
        t_max = uta.trajectory[-1].time
    if step < 1:
        step = 1
    dt = np.round((uta.trajectory[1].time - uta.trajectory[0].time),3)
    frame_ids = np.arange(int((t_min - uta.trajectory[0].time)/dt), int((t_max - uta.trajectory[0].time)/dt + 1), step)
    dt = dt*step

    r1 = []; cells = []
    for frame in frame_ids:

        if frame%5000 == 0:
            print("Frame " + str(frame))
        
        ts = uta.trajectory[frame]
        cell = ts.dimensions
        cells.append(cell)

        r1.append(a1.positions)

    with h5py.File('/tmp/Cluster.hdf5','w') as f:
        dset1 = f.create_dataset("cells", data = cells)
        dset2 = f.create_dataset("r1", data=r1)
        dset4 = f.create_dataset("frames", data=frame_ids)



def main(trj_file, top_file, aw_name, coord, t_min, t_max, step, nt):

    # Load in the trajectory file
    if not os.path.exists('/tmp/Cluster.hdf5'):
        load_TRR()
        exit()

    with h5py.File('/tmp/Cluster.hdf5','r') as f:
        dset1 = f['frames']; frame_ids = dset1[:]

    print("Pairing Analysis")
    pool = mp.Pool(processes=nt)
    func = functools.partial(iPair)
    return_arr = pool.map(func, range(len(frame_ids)))
    pool.close()
    pool.join()
    return_arr = np.array(return_arr)

    clust = np.mean(return_arr[:,0,:], axis = 0)
    while clust[-1] == 0:
        clust = clust[:-1]
    
    max_dist = np.sum(return_arr[:,1,:], axis = 0); count = np.sum(return_arr[:,2,:], axis = 0)
    max_dist = np.divide(max_dist, count, out=np.zeros_like(max_dist), where=count!=0)
    while max_dist[-1] == 0:
        max_dist = max_dist[:-1]
    
    S_max_std = np.std(return_arr[:,3,0], axis = 0)
    S_max = np.mean(return_arr[:,3,0], axis = 0)
    with open('S_max.dat', 'w') as anaout:
        print('# S_max', file=anaout)
        print(' 0.00 {:10.5f} {:1.5e}'.format(S_max, S_max_std), file=anaout)
    
    P_frac_std = np.std(return_arr[:,4,0])
    P_frac = np.mean(return_arr[:,4,0])
    with open('P_frac.dat', 'w') as anaout:
        print('# Percolation Fraction', file=anaout)
        print(' 0.00 {:10.5f} {:1.5e}'.format(P_frac, P_frac_std), file=anaout)

    pairPrint2('Clust_water', clust)
    pairPrint3('Size_water', max_dist)
    
    # Deletes the temporary h5py file
    os.remove('/tmp/Cluster.hdf5')

if __name__ == "__main__":
    main(trj_file, top_file, aw_name, coord, t_min, t_max, step, nt)
