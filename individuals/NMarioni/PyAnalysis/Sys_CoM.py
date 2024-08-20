#This script saves the System Center of Mass (CoM) to a file. This is so you can still retrieve the correct CoM using reduced .trr/.xtc files that DONT contain all atoms

import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
import numpy as np
import h5py# h5py is used for efficient array storage and to free up memory by using it for I/O

import multiprocessing as mp
import functools
import os

from sys import argv
script, trj_file, top_file = argv
# trj_file = unwrapped .trr/.xtc file
#   NOTE: Must contain ALL atoms and be unwrapped
# top_file = .tpr file

#Example: python3 ${Path}/Sys_CoM.py unwrap.trr md.tpr

def load_TRR():
# Load in the trajectory file and write necessary data to a h5py file
# For very large files, it is recommended to dump only the atoms you are interested in to a separate .trr/.xtc file for analysis

    global t_min, t_max, step

    uta = mda.Universe(top_file, trj_file)

    t_min = uta.trajectory[0].time
    t_max = uta.trajectory[-1].time
    step = 1
    dt = np.round((uta.trajectory[1].time - uta.trajectory[0].time),3)
    frame_ids = np.arange(int((t_min - uta.trajectory[0].time)/dt), int((t_max - uta.trajectory[0].time)/dt + 1), step)
    dt = dt*step
    print("Timestep " + str(dt))
    
    times = []; CoM = []
    for frame in frame_ids:
        ts = uta.trajectory[frame]
    
        if ts.time%5000 == 0:
            print("Time "+ str(ts.time))
        
        times.append(ts.time)
        CoM.append(uta.atoms.center_of_mass())

    times = np.array(times); CoM = np.array(CoM)

    with h5py.File('CoM.hdf5','w') as f:
        dset1 = f.create_dataset("times", data=times)
        dset2 = f.create_dataset("CoM", data=CoM)



def main(trj_file, top_file):
    load_TRR()



if __name__ == "__main__":
    main(trj_file, top_file)
