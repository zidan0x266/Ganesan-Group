# This code generates a hydrogen-bonding probability distribution

import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
import MDAnalysis.lib.distances as distances
import numpy as np

import multiprocessing as mp
import functools
import h5py
import os
from sys import argv

script, trj_file, top_file, a1_name, cutoff, t_min, t_max, step, nt = argv
cutoff = float(cutoff); t_min = float(t_min); t_max = float(t_max); step = int(step); nt = int(nt)
# trj_file = .trr/.xtc file, top_file = .tpr file
# a1_name = name of atom hydrogen bonding to
# cutoff = coordination distance between a1 and a2 (Angstroms); 3.5 angstroms by default
# t_min = start time for analysis in ps (-1 assumes the start time of the first frame)
# t_max = end time for analysis in ps (-1 assumes the end time of the last frame)
# step = frame step-size (-1 assumes a step-size of 1)
# nt = number of threads

if cutoff == -1:
    cutoff = 3.5



def pairPrint(filename, pair):# print ion pair proability distribution
    with open('{}.dat'.format(filename), 'w') as anaout:
        print('# HBonds Probability', file=anaout)
        if a1_name == "OW":
            print('# Number of H-Bonds per water molecule (includes HBonds to O and from H', file=anaout)
        for i in range(0, len(pair)):
            print('{:10.3f} {:10.5f}'.format(i, pair[i]), file=anaout)



def iPair(frame):
# Determine the hydrogen bonds formed within this frame

    if frame%5000 == 0:
        print("Frame " + str(frame))

    with h5py.File('/tmp/HBond.hdf5','r') as f:
        dset1 = f['r1']; a1 = dset1[frame] # O Acceptor you are analyzing
        dset2 = f['r2']; a2 = dset2[frame] # Water Ow
        dset3 = f['r3']; a3 = dset3[frame] # Water Hw1
        dset4 = f['r4']; a4 = dset4[frame] # Water Hw2
        dset3 = f['cells']; cell = dset3[frame]

    # Check O - Ow distance less than cutoff (3.5 Angstroms by default)
    distpair = distances.capped_distance(a1, a2, cutoff, box=cell)[0]

    hbond = np.zeros(len(a1))
    for pair in distpair:
        #Skip self interactions for water-water HBonding (if O Acceptor = Ow)
        if a1_name == "OW" and pair[0] == pair[1]:
            continue

        #Distance between target and potential HBond donors
        h_dist = distances.distance_array(a1[pair[0]], np.array([a3[pair[1]], a4[pair[1]]]), box=cell)[0]
        
        ba = a1[pair[0]] - a2[pair[1]]
        # Guess HBond arises due to the closest H
        if h_dist[0] < h_dist[1]:
            bc = a3[pair[1]] - a2[pair[1]]
            bc_temp = a4[pair[1]] - a2[pair[1]]
        else:
            bc = a4[pair[1]] - a2[pair[1]]
            bc_temp = a3[pair[1]] - a2[pair[1]]
        
        angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        # Deal with math errors
        if angle > 1 and angle < 1.0001:
            angle = np.arccos(1) * (180/np.pi)
        elif angle < -1 and angle > -1.0001:
            angle = np.arccos(-1) * (180/np.pi)
        else:
            angle = np.arccos(angle) * (180/np.pi)

        ##If np.accros raises an error
        #if np.isnan(angle):
        #    print("Is NaN")
        #    print(np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc)))

        if angle < 30:
            hbond[pair[0]] += 1

            #Number of HBonds to water includes accepted and donated HBonds
            if a1_name == "OW":
                hbond[pair[1]] += 1

        # If closest H not HBond, check other H
        else:
            bc = bc_temp

            angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
            # Deal with math errors
            if angle > 1 and angle < 1.0001:
                angle = np.arccos(1) * (180/np.pi)
            elif angle < -1 and angle > -1.0001:
                angle = np.arccos(-1) * (180/np.pi)
            else:
                angle = np.arccos(angle) * (180/np.pi)

            ##If np.accros raises an error
            #if np.isnan(angle):
            #    print("Is NaN")
            #    print(np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc)))

            if angle < 30:
                hbond[pair[0]] += 1

                #Number of HBonds to water includes accepted and donated HBonds
                if a1_name == "OW":
                    hbond[pair[1]] += 1

    return hbond



def pair_dist(pair_ar):
# Turn the array of hydrogen bonds into a simple 1-D distribution

    with h5py.File('/tmp/HBond.hdf5','r') as f:
        dset1 = f['r1']; a1 = len(dset1[0])

    npair = np.zeros(10)# Size of the distribution. Can adjust as needed
    for pair in pair_ar:
        for i in range(len(npair)):
            npair[int(i)] += np.count_nonzero(pair == i)

    return npair / a1 / len(pair_ar)



def load_TRR():
# loads in the trajectory and saves the necessary data to a temporary h5py file

    global t_min, t_max, step

    uta = mda.Universe(top_file, trj_file)

    a1 = uta.select_atoms("name " + a1_name)
    a2 = uta.select_atoms("name OW")
    a3 = uta.select_atoms("name HW1")
    a4 = uta.select_atoms("name HW2")

    if t_min == -1:
        t_min = uta.trajectory[0].time
    if t_max == -1:
        t_max = uta.trajectory[-1].time
    if step < 1:
        step = 1
    dt = np.round((uta.trajectory[1].time - uta.trajectory[0].time),3)
    frame_ids = np.arange(int((t_min - uta.trajectory[0].time)/dt), int((t_max - uta.trajectory[0].time)/dt + 1), step)
    dt = dt*step
    print("Timestep " + str(dt))

    r1 = []; r2 = []; r3 = []; r4 = []; cells = []
    for frame in frame_ids:

        if frame%5000 == 0:
            print("Frame " + str(frame))

        ts = uta.trajectory[frame]
        cell = ts.dimensions
        cells.append(cell)

        r1.append(a1.positions)
        r2.append(a2.positions)
        r3.append(a3.positions)
        r4.append(a4.positions)

    with h5py.File('/tmp/HBond.hdf5','w') as f:
        dset1 = f.create_dataset("r1", data=r1)
        dset2 = f.create_dataset("r2", data=r2)
        dset3 = f.create_dataset("r3", data=r3)
        dset4 = f.create_dataset("r4", data=r4)
        dset5 = f.create_dataset("cells", data = cells)
        dset6 = f.create_dataset("dt", data = [dt])
        dset7 = f.create_dataset("frames", data = frame_ids)



def main(trj_file, top_file, a1_name, cutoff, t_min, t_max, step, nt):

    # Makes sure the temporary h5py file does not exist
    #if os.path.exists('/tmp/HBond.hdf5'):
    #    os.remove('/tmp/HBond.hdf5')

    # Load in the trajectory file
    if not os.path.exists('/tmp/HBond.hdf5'):
        load_TRR()
        exit()

    with h5py.File('/tmp/HBond.hdf5','r') as f:
        dset1 = f['frames']; frame_ids = dset1[:]
    frame_ids = frame_ids - frame_ids[0]

    # Determine the number of bound a2 atoms for each a1 atom at each frame
    print("Pairing Analysis")
    pool = mp.Pool(processes=nt)
    func = functools.partial(iPair)
    pair_temp = pool.map(func, range(len(frame_ids)))
    pool.close()
    pool.join()

    # Probability of an a1 atom being bound to x a2 atoms
    print("Probability Dist")
    npair = pair_dist(pair_temp)
    pairPrint("HBond_{}".format(a1_name), npair)

    # Deletes the temporary h5py file
    os.remove('/tmp/HBond.hdf5')

if __name__ == "__main__":
    main(trj_file, top_file, a1_name, cutoff, t_min, t_max, step, nt)
