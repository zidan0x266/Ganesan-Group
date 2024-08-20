#This script calculates the 1-D Displacement over time (DIS)

import MDAnalysis as mda
import MDAnalysis.lib.distances as distances
import numpy as np
import h5py

import multiprocessing as mp
import functools
import os

from sys import argv
script, trj_file, top_file, a_name, dimension, t_min, t_max, step, nt, CoM_check = argv
t_min = float(t_min); t_max = float(t_max); step = int(step); nt = int(nt)
# trj_file = .trr/.xtc file, top_file = .tpr file
# a_name = name of atom to be analyzed
# dimension = string containing 'x', 'y', or 'z', indicating the axis to calculate displacement on
# t_min = start time for analysis in ps (-1 assumes the start time of the first frame)
# t_max = end time for analysis in ps (-1 assumes the end time of the last frame)
# step = frame step-size (-1 assumes a step-size of 1)
# nt = number of threads
# NOTE: Make sure to have unwrapped coordinates (gmx trjconv with -pbc nojump)
# CoM_check = 0 if the trj_file contains ALL atoms, 1 else
# NOTE: If the trj_file does NOT contain ALL atoms, then the system center of mass will not be correct. First run Sys_CoM.py and then send CoM_check = 1

# Example: python3 ${Path}/analysis_dis.py unwrap.trr md.tpr NA x -1 -1 -1 128 0

dim_ar = np.array(['x', 'y', 'z'])
dimension = np.where(dim_ar == dimension)[0][0]



def dis_analysis(dim, dt, dt_max):
# calculate the DIS with CoM removed
#
# Inputs: dim = index of the desired dimension, dt = timestep in ps; dt_max = max time to calculate DIS over

    nframe = int(round(dt_max/dt))

    pool = mp.Pool(processes=nt)
    func = functools.partial(dis_calc, dim, nframe)
    dis = pool.map(func, list(range(nframe)))
    pool.close()
    pool.join()

    dis_bins = np.arange(0,nframe*dt,dt)
    dis = np.array(dis)

    with open('dis_{}.xvg'.format(a_name), 'w') as anaout:
        print('# Time DIS (nm)', file=anaout)
        for i, bin_i in enumerate(dis_bins):
            if i >= len(dis):
                print('{:10.3f} {:10.5f}'.format(bin_i, 0.0), file=anaout)
            else:
                print('{:10.3f} {:10.5f}'.format(bin_i, dis[i]), file=anaout)



def dis_calc(dim, nframe, df):
# calculate the DIS with CoM removed
#
# Inputs: dim = index of the desired dimension, nframe = total number of frames, df = frame step to calculate the DIS over

    if df == 0:
        return 0.0

    if df%5000 == 0:
        print("dFrame "+ str(df))

    DIS_file = h5py.File('/tmp/r_DIS_'+a_name+'.hdf5','r'); r = DIS_file['r'][:,:,dim]

    dis = 0.0; count = 0
    for j in range(0, nframe-1, 5):# Statistical enchancement by calculating over every 10 frames
        if df + j >= nframe:
            break
        
        dis += np.mean(r[j+df,:] - r[j,:])
        count += 1

    DIS_file.close()

    return dis / count



def load_TRR():
# Load in the trajectory file and write necessary data to a h5py file
# For very large files, it is recommended to dump only the atoms you are interested in to a separate .trr/.xtc file for analysis
# NOTE: Must supply with unwrapped coordinates

    global t_min, t_max, step

    uta = mda.Universe(top_file, trj_file, tpr_resid_from_one=True)
    a = uta.select_atoms("name " + a_name)

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

    if CoM_check == '1':
        with h5py.File('CoM.hdf5','r') as f:
            times = f['times'][:]
            CoM_ar = f['CoM'][:,:]
    
    r = []
    for frame in frame_ids:
        ts = uta.trajectory[frame]
        cell = ts.dimensions

        if CoM_check == '0':
            CoM = uta.atoms.center_of_mass()
        elif CoM_check == '2':
            CoM = uta.select_atoms("moltype MOL").center_of_mass()
        else:
            CoM = CoM_ar[np.where(times == ts.time)]
    
        if ts.time%5000 == 0:
            print("Time "+ str(ts.time))

        r.append((a.positions - CoM)/10)

    with h5py.File('/tmp/r_DIS_'+a_name+'.hdf5','w') as f:
        dset1 = f.create_dataset("r", data=r)
        dset2 = f.create_dataset("dt", data=[dt])
        dset3 = f.create_dataset("frames", data=frame_ids)



def main(trj_file, top_file, a_name, dimension, t_min, t_max, step, nt, CoM_check):

    # If there is not a position h5py file, then create one and end the program
    # This is done to avoid memory problems during multiprocessing
    if not os.path.exists('/tmp/r_DIS_'+a_name+'.hdf5'):
        load_TRR()
        exit()
    with h5py.File('/tmp/r_DIS_'+a_name+'.hdf5','r') as f:
        dset1 = f['dt']; dt = dset1[0]
        dset2 = f['frames']; frame_ids = dset2[:]
   
    print("DIS Analysis")
    dt_max = len(frame_ids)*dt# as written dt_max is the total time frame
    #dt_max=10000
    dis_analysis(dimension, dt, dt_max)

    os.remove('/tmp/r_DIS_'+a_name+'.hdf5')



if __name__ == "__main__":
    main(trj_file, top_file, a_name, dimension, t_min, t_max, step, nt, CoM_check)
