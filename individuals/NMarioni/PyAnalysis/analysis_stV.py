# Original code written by Everett Zofchak, adapted to vehicular analysis by Nico Marioni
# Calculates the ion vehicular autcorrelation funtions (continuous, stV, and/or dicontinuous, ctV)

import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
import MDAnalysis.lib.distances as distances
import numpy as np

import multiprocessing as mp
import functools
import h5py
import os
from sys import argv

#from memory_profiler import profile
## Memory_profiler is useful for dealing with memory allocation errors
## To use, add the "@profile" decorator to the line above the function you want to see the line-by-line memory usage of
##   NOTE: The code must run to completion for the memory_profiler to print without error. Stopping the code with Ctrl+C or exit() will not work with it


script, trj_file, top_file, a1_name, a2_name, coord, t_min, t_max, step, t_calc, nt = argv
coord = float(coord); t_min = float(t_min); t_max = float(t_max); step = int(step); t_calc = float(t_calc); nt = int(nt)
# trj_file = .trr/.xtc file, top_file = .tpr file
# a1_name = name of reference atom to be analyzed
# a2_name = name of selection atom to be analyzed
# coord = coordination distance between a1 and a2 (Angstroms)
# t_min = start time for analysis in ps (-1 assumes the start time of the first frame)
# t_max = end time for analysis in ps (-1 assumes the end time of the last frame)
# step = frame step-size (-1 assumes a step-size of 1)
# t_calc = time in ps to calculate c_phiphi out to (-1 assumes t_min - t_max)
# nt = number of threads

# Example: python3 ${Path}/analysis_stV.py md.trr md.tpr NA CL 3.67 75000 100000 1 10000 128



def stct_calc1(num_s, num_c, num_count, df_max, df):
# Calculates both the continuous and discontinuous vehicular autocorrelation function
#
# Inputs: num_s = array containing st vehicular data; num_c = array containing ct vehicular data; num_count = array counting frames analyzed
#         df_max = max frame step for calculations; df = current frame step
        if df%1000 == 0:
                print("dFrame "+ str(df))

        f = h5py.File('/tmp/stct.hdf5','r'); dset2 = f['master']

        arr_s = dset2[df]; save = dset2[df]; Veh = np.sum(save, axis = 1)
        num_s[0] += np.sum(arr_s); num_c[0] += np.sum(save); num_count[0] += 1
        for i in range(df+1, min(df+1+df_max, len(dset2))):
                arr = dset2[i]

                arr_s = np.multiply(arr_s, arr); arr_c = np.multiply(save, arr)
                num_s[i-df] += np.sum(np.divide(np.sum(arr_s, axis=1), Veh, where=(Veh!=0), out=np.zeros_like(Veh,dtype=np.float32)))
                num_c[i-df] += np.sum(np.divide(np.sum(arr_c, axis=1), Veh, where=(Veh!=0), out=np.zeros_like(Veh,dtype=np.float32)))
                num_count[i-df] += 1

        f.close()

        return num_s, num_c, num_count

def stct_calc2(num_c, num_count, df_max, df):
# Calculates the discontinuous vehicular autocorrelation function
#
# Inputs: num_c = array containing ct vehicular data; num_count = array counting frames analyzed
#         df_max = max frame step for calculations; df = current frame step
        if df%1000 == 0:
                print("dFrame "+ str(df))

        f = h5py.File('/tmp/stct.hdf5','r'); dset2 = f['master']

        save = dset2[df]; Veh = np.sum(save, axis = 1)
        num_c[0] += np.sum(save); num_count[0] += 1
        for i in range(df+1, min(df+1+df_max, len(dset2))):
                arr = dset2[i]

                arr_c = np.multiply(save, arr)
                num_c[i-df] += np.sum(np.divide(np.sum(arr_c, axis=1), Veh, where=(Veh!=0), out=np.zeros_like(Veh,dtype=np.float32)))
                num_count[i-df] += 1

        f.close()

        return num_c, num_count

def stct_calc3(num_s, num_count, df_max, df):
# Calculates the continuous vehicular autocorrelation function
#
# Inputs: num_s = array containing st vehicular data; num_count = array counting frames analyzed
#         df_max = max frame step for calculations; df = current frame step
        if df%1000 == 0:
                print("dFrame "+ str(df))

        f = h5py.File('/tmp/stct.hdf5','r'); dset2 = f['master']

        arr_s = dset2[df]; save = dset2[df]; Veh = np.sum(save, axis = 1)
        num_s[0] += np.sum(arr_s); num_count[0] += 1
        for i in range(df+1, min(df+1+df_max, len(dset2))):
                arr = dset2[i]

                arr_s = np.multiply(arr_s, arr)
                num_s[i-df] += np.sum(np.divide(np.sum(arr_s, axis=1), Veh, where=(Veh!=0), out=np.zeros_like(Veh,dtype=np.float32)))
                num_count[i-df] += 1

        f.close()

        return num_s, num_count



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
    if step == -1:
        step = 1
    
    dt = np.round((uta.trajectory[1].time - uta.trajectory[0].time),3)
    frame_ids = np.arange(int((t_min - uta.trajectory[0].time)/dt), int((t_max - uta.trajectory[0].time)/dt + 1), step)
    dt = dt*step
    print("Timestep " + str(dt))

    master = []
    for frame in frame_ids:

        if frame%1000 == 0:
            print("Frame " + str(frame))
    
        ts = uta.trajectory[frame]
        cell = ts.dimensions

        r1 = a1.positions
        r2 = a2.positions

        dists = distances.capped_distance(r1,r2,coord,box=cell)[0]
        coord_mat = np.zeros((len(r1),len(r2)),dtype=bool)
        for i in dists:
            coord_mat[i[0],i[1]] = True
        master.append(coord_mat)
    
    #master = np.array(master)

    with h5py.File('stct.hdf5','w') as f:
        dset1 = f.create_dataset("dt", data = [dt])
        dset2 = f.create_dataset("frames", data = frame_ids)
        dset3 = f.create_dataset("master", data = master)



def main(trj_file, top_file, a1_name, a2_name, coord, t_min, t_max, step, t_calc, nt):

    # Makes sure the temporary h5py file does not exist

    # Load in the trajectory file
    if not os.path.exists('stct.hdf5'):
        load_TRR()
        exit()

    # Checks file size and copies to /tmp/
    if not os.path.exists('/tmp/stct.hdf5'):
        file_size = os.path.getsize('stct.hdf5') / (1024**3)
        print(file_size)
        if file_size > 125.0:
            print("File too large: " + file_size + "GB")
            exit()
        else:
            os.system('cp stct.hdf5 /tmp/')
    
    with h5py.File('/tmp/stct.hdf5','r') as f:
        dset1 = f['frames']; frame_ids = dset1[:]
        dset2 = f['dt']; dt = dset2[0]

    if t_calc == -1:
        if step == -1:
            dt_max = -(frame_ids[-1] - frame_ids[0])*dt/step # Time to calculate st out to in ps (10000 is 10 ns)
        else:
            dt_max = (frame_ids[-1] - frame_ids[0])*dt/step # Time to calculate st out to in ps (10000 is 10 ns)
    else:
        dt_max = t_calc
    print(dt_max)





    # Calculates both the continuous and discontinuous vehicular autocorrelation function
    num_s, num_c, num_count = np.zeros(int(dt_max / dt)+1,dtype=np.float32), np.zeros(int(dt_max / dt)+1,dtype=np.float32), np.zeros(int(dt_max / dt)+1,dtype=np.uintc)
    pool = mp.Pool(processes=nt)
    resevoir = [ti for ti in range(0,len(frame_ids),5)]
    func = functools.partial(stct_calc1,num_s,num_c,num_count, int(dt_max / dt))
    return_data = pool.map(func, resevoir)
    pool.close()
    pool.join()
    return_data = np.array(return_data)
    
    num_s = np.sum(return_data[:,0,:], axis = 0)
    num_c = np.sum(return_data[:,1,:], axis = 0)
    num_count = np.sum(return_data[:,2,:], axis = 0)
    
    st = np.divide(num_s, num_count, where=num_count!=0); st = np.divide(st, st[0], where=st[0]!=0)
    ct = np.divide(num_c, num_count, where=num_count!=0); ct = np.divide(ct, ct[0], where=ct[0]!=0)

    print("Printing Data to Files")
    with open('stV_{}_{}.xvg'.format(a1_name, a2_name), 'w') as anaout:
        print('# Time st(t)', file=anaout)
        for i in range(len(st)):
            print('{:10.3f} {:10.5f}'.format(i*dt, st[i]), file=anaout)

    with open('ctV_{}_{}.xvg'.format(a1_name, a2_name), 'w') as anaout:
        print('# Time ct(t)', file=anaout)
        for i in range(len(ct)):
            print('{:10.3f} {:10.5f}'.format(i*dt, ct[i]), file=anaout)

    os.remove('stct.hdf5'); os.remove('/tmp/stct.hdf5')





    ## Calculates the discontinuous vehicular autocorrelation function
    #num_c, num_count = np.zeros(int(dt_max / dt)+1,dtype=np.float32), np.zeros(int(dt_max / dt)+1,dtype=np.uintc)
    #pool = mp.Pool(processes=nt)
    #resevoir = [ti for ti in range(0,len(frame_ids),5)]
    #func = functools.partial(stct_calc2,num_c,num_count, int(dt_max / dt))
    #return_data = pool.map(func, resevoir)
    #pool.close()
    #pool.join()
    #return_data = np.array(return_data)

    #num_c = np.sum(return_data[:,0,:], axis = 0)
    #num_count = np.sum(return_data[:,1,:], axis = 0)

    #ct = np.divide(num_c, num_count, where=num_count!=0); ct = np.divide(ct, ct[0], where=ct[0]!=0)

    #print("Printing Data to Files")
    #with open('ctV_{}_{}.xvg'.format(a1_name, a2_name), 'w') as anaout:
    #    print('# Time ct(t)', file=anaout)
    #    for i in range(len(ct)):
    #        print('{:10.3f} {:10.5f}'.format(i*dt, ct[i]), file=anaout)

    #os.remove('stct.hdf5'); os.remove('/tmp/stct.hdf5')





    ## Calculates the continuous vehicular autocorrelation function
    #num_s, num_count = np.zeros(int(dt_max / dt)+1,dtype=np.float32), np.zeros(int(dt_max / dt)+1,dtype=np.uintc)
    #pool = mp.Pool(processes=nt)
    #resevoir = [ti for ti in range(0,len(frame_ids),5)]
    #func = functools.partial(stct_calc3,num_s,num_count, int(dt_max / dt))
    #return_data = pool.map(func, resevoir)
    #pool.close()
    #pool.join()
    #return_data = np.array(return_data)
    #
    #num_s = np.sum(return_data[:,0,:], axis = 0)
    #num_count = np.sum(return_data[:,1,:], axis = 0)
    #
    #st = np.divide(num_s, num_count, where=num_count!=0); st = np.divide(st, st[0], where=st[0]!=0)

    #print("Printing Data to Files")
    #with open('stV_{}_{}.xvg'.format(a1_name, a2_name), 'w') as anaout:
    #    print('# Time st(t)', file=anaout)
    #    for i in range(len(st)):
    #        print('{:10.3f} {:10.5f}'.format(i*dt, st[i]), file=anaout)

    #os.remove('stct.hdf5'); os.remove('/tmp/stct.hdf5')
    
if __name__ == "__main__":
    main(trj_file, top_file, a1_name, a2_name, coord, t_min, t_max, step, t_calc, nt)
