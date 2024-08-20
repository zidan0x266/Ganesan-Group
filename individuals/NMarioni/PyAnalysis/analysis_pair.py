# This script calculates the ion pairing probability distribution

import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
import MDAnalysis.lib.distances as distances
import numpy as np

import multiprocessing as mp
import functools
import h5py
import os
from sys import argv

script, trj_file, top_file, a1_name, a2_name, coord, t_min, t_max, step, nt = argv
coord = float(coord); t_min = float(t_min); t_max = float(t_max); step = int(step); nt = int(nt)
# trj_file = .trr/.xtc file, top_file = .tpr file
# a1_name = name of reference atom to be analyzed
# a2_name = name of selection atom to be analyzed
# coord = coordination distance between a1 and a2
# t_min = start time for analysis in ps (-1 assumes the start time of the first frame)
# t_max = end time for analysis in ps (-1 assumes the end time of the last frame)
# step = frame step-size (-1 assumes a step-size of 1)
# nt = number of threads



def pairPrint(filename, pair):  # print displacement
    with open('{}.dat'.format(filename), 'w') as anaout:
        print('# ' + filename + ' Probability', file=anaout)
        for i in range(0, len(pair)):
            print('{:10.3f} {:10.5f}'.format(i, pair[i]), file=anaout)



def iPair(frame):
# finds number of a2 per a1 at the given frame
#
# Inputs: frame = time frame to be analyzed

    with h5py.File('/tmp/Pair.hdf5','r') as f:
        dset1 = f['r1']; a1 = dset1[frame]
        dset2 = f['r2']; a2 = dset2[frame]
        dset3 = f['cells']; cell = dset3[frame]

    if frame%500 == 0:
        print("Frame " + str(frame))

    pair = np.zeros(10)
    for a1_i in a1:
        a2_tmp = a2
        
        if a1_name == a2_name:
            a2_tmp = np.delete(a2_tmp, np.where((a2_tmp == a1_i).all(axis=1))[0][0], 0)
        
        distpair = distances.capped_distance(a1_i, a2_tmp, coord, box=cell)[1]

        pair[len(distpair)] += 1

    return pair / len(a1)



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
    print("Timestep " + str(dt))

    r1 = []; r2 = []; cells = []
    for frame in frame_ids:
        ts = uta.trajectory[frame]
        cells.append(ts.dimensions)

        r1.append(a1.positions)
        r2.append(a2.positions)

    cells = np.array(cells); r1 = np.array(r1); r2 = np.array(r2)

    with h5py.File('/tmp/Pair.hdf5','w') as f:
        dset1 = f.create_dataset("r1", data=r1)
        dset2 = f.create_dataset("r2", data=r2)
        dset3 = f.create_dataset("cells", data = cells)

    return frame_ids



def main(trj_file, top_file, a1_name, a2_name, coord, t_min, t_max, step, nt):
    frame_ids = load_TRR()

    pool = mp.Pool(processes=nt)
    func = functools.partial(iPair)
    pair_ar = pool.map(func, list(frame_ids))
    pool.close()
    pool.join()

    pairPrint("aPair_{}_{}".format(a1_name, a2_name), np.mean(pair_ar, axis = 0))

    os.remove('/tmp/Pair.hdf5')

if __name__ == "__main__":
    main(trj_file, top_file, a1_name, a2_name, coord, t_min, t_max, step, nt)
