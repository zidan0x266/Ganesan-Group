# Calculates the static and dynamic structure factor

import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
import MDAnalysis.lib.distances as distances
import numpy as np
import random

import multiprocessing as mp
import functools
import h5py
import os
from sys import argv


script, trj_file, top_file, poly_names, t_min, t_max, step, t_calc, nt = argv
t_min = float(t_min); t_max = float(t_max); step = int(step); t_calc = float(t_calc); nt = int(nt)
# trj_file = .trr/.xtc file, top_file = .tpr file
# poly_names = string of polymer residues, e.g. resname SPS or resname DVB
# t_min = start time for analysis in ps (-1 assumes the start time of the first frame)
# t_max = end time for analysis in ps (-1 assumes the end time of the last frame)
# step = frame step-size (-1 assumes a step-size of 1)
# t_calc = time in ps to calculate sqt out to (-1 assumes t_min - t_max)
# nt = number of threads

# Example: python3 ${Path}/analysis_sq.py md.trr md.tpr 'resname PEO or resname ...' -1 -1 -1 5000 128



def sqPrint(filename, sq):  # print sq to the output
    with open('{}.xvg'.format(filename), 'w') as sqtout:
        print('# q S(q)', file=sqtout)
        for i in range(len(sq[0])):
            print('{:10.3f} {:10.5f}'.format(sq[0][i], sq[1][i]), file=sqtout)



def sqtPrint(filename, sqt, dt):  # print sqt to the output
    with open('{}.xvg'.format(filename), 'w') as sqtout:
        print('# time 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0', file=sqtout)
        for i in range(len(sqt[0])):
            print('{:10.3f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(i*dt, 
            sqt[0][i], 
            sqt[1][i], 
            sqt[2][i], 
            sqt[3][i], 
            sqt[4][i], 
            sqt[5][i], 
            sqt[6][i], 
            sqt[7][i], 
            sqt[8][i], 
            sqt[9][i], 
            sqt[10][i],
            sqt[11][i],
            sqt[12][i],
            sqt[13][i],
            sqt[14][i],
            sqt[15][i],
            sqt[16][i],
            sqt[17][i],
            sqt[18][i],
            sqt[19][i],
            sqt[20][i],), file=sqtout)



def presq(cell):
    """
    This function generates the number of q values, since the number
    depends on the box length (cell), the structure is similar to 
    the following qgenerator function
    """
    minq = 2 * np.pi / cell # 2pi/L
    largest_q = 5.05 # calculates sq out to q = 5 Angstroms^-1
    maximum_q = 100 # 100 q vectors per bin of q
    n_max = int(np.ceil(largest_q/minq))

    num_qvec = np.zeros(1000); counter = 0
    for nx in range(-n_max, n_max+1, 1):
        for ny in range(-n_max, n_max+1, 1):
            for nz in range(-n_max, n_max+1, 1):

                if nx == 0 and ny == 0 and nz == 0:
                    continue

                qx = nx * minq; qy = ny * minq; qz = nz * minq
                q_sqrt = np.sqrt(qx**2 + qy**2 + qz**2)

                if q_sqrt > largest_q:
                    continue

                remain = (int(q_sqrt * 1000)) % 100
                if remain >= 90 or remain <= 10:
                    qvalue_10 = int(round(q_sqrt * 10))  # q values multiply by 10

                    if num_qvec[qvalue_10] < maximum_q:                    
                        counter += 1
                        num_qvec[qvalue_10] += 1

    return counter



def sqgenerator(number_q, cell):
    """
    This function generates the q vector
    """
    total_q = number_q
    minq = 2 * np.pi / cell  # 2pi/L
    largest_q = 5.05 # calculates sq out to q = 5 Angstroms^-1
    maximum_q = 100 # 100 q vectors per bin of q
    n_max = int(np.ceil(largest_q/minq))

    q = np.zeros(total_q); qtmps = []
    q_vec = np.zeros((total_q, 3)); q_vec_tmps = []
    for nx in range(-n_max, n_max+1, 1):
        for ny in range(-n_max, n_max+1, 1):
            for nz in range(-n_max, n_max+1, 1):

                if nx == 0 and ny == 0 and nz == 0:
                    continue

                qx = nx * minq; qy = ny * minq; qz = nz * minq
                q_sqrt = np.sqrt(qx**2 + qy**2 + qz**2)

                if q_sqrt > largest_q:
                    continue

                remain = (int(q_sqrt * 1000)) % 100  # assure the generated q has a variance of +- 10%
                if remain >= 90 or remain <= 10:
                    qtmps.append(q_sqrt)
                    q_vec_tmps.append([qx,qy,qz])

    # randomly selects maximum_q number of q_vectors for each bin of q
    qtmps = np.array(qtmps); q_vec_tmps = np.array(q_vec_tmps); count = 0
    for i in range(51):
        arr = np.where((qtmps > 0.1*i - 0.05) & (qtmps < 0.1*i + 0.05))[0]

        if len(arr) > maximum_q:
            arr = np.array(random.sample(list(arr), k=maximum_q))

        for j in arr:
            q[count] = qtmps[j]
            q_vec[count] = q_vec_tmps[j]
            count += 1

    return q_vec, q, total_q



def calc_sq(nframes, scatter_length, qs):
# Calculates sq based on q_vec, qs
# Inputs: nframes = total number of frames, scatter_length = scatter factor of each atom; qs = q_vector
    f = h5py.File('/tmp/sq.hdf5','r'); dset1 = f['master']

    sq_arr = np.zeros(nframes)
    for i in range(nframes):
        pos = dset1[i]

        qr = (qs * pos).sum(axis=1)
        cossum = (scatter_length * np.cos(qr)).sum()
        sinsum = (scatter_length * np.sin(qr)).sum()

        sq_arr[i] += cossum **2 + sinsum **2
    
    f.close()

    return sq_arr



def sq_analysis(cell, nframes, scatter_length):
# Calculate the static structure factor
# inputs: cell = simulation cell length; nframes = total number of frames; scatter_length = scatter factor of each atom

    # Get q vectors (q = 2pi/L*(nx, ny, nz), where ni is an integrer)
    number_q = presq(cell)
    q_vec, q, n_max = sqgenerator(number_q, cell)

    q_points = 51 # determined by largest_q, i.e., (5.05 + 0.05)*10
    sqtmps = np.zeros((q_points, nframes)); count = np.zeros((q_points))
    sqvals = np.zeros((2, q_points))  # bins from 0.0 to 5.0 with increament 0.1

    # calculate sq at each q_vector
    pool = mp.Pool(processes=nt)
    fn_q_analysis = functools.partial(calc_sq, nframes, scatter_length)
    sq = pool.map(fn_q_analysis, q_vec); sq = np.array(sq)
    pool.close()
    pool.join()

    # bin sq based on the value of q
    for qindex in range(n_max):
        sqtmps[int(round(q[qindex] * 10))] += sq[qindex]
        count[int(round(q[qindex] * 10))] += nframes
    sq_sum = np.sum(sqtmps / np.square(scatter_length).sum(), axis = 1)
    
    # normallization
    for i in range(1, q_points):
        sqvals[0][i] = i * 0.1
        if count[i] != 0:
            sqvals[1][i] = sq_sum[i] / count[i]

    sqPrint('sq', sqvals)



def presqt(cell):
    """
    This function generates the number of q values, since the number
    depends on the box length (cell), the structure is similar to 
    the following qgenerator function
    """
    minq = 2 * np.pi / cell # 2pi/L
    largest_q = 2.05 # calculates sqt out to q = 2 Angstroms^-1
    maximum_q = 100 # 100 q vectors per bin of q
    n_max = int(np.ceil(largest_q/minq))

    num_qvec = np.zeros(1000); counter = 0
    for nx in range(-n_max, n_max+1, 1):
        for ny in range(-n_max, n_max+1, 1):
            for nz in range(-n_max, n_max+1, 1):

                if nx == 0 and ny == 0 and nz == 0:
                    continue

                qx = nx * minq; qy = ny * minq; qz = nz * minq
                q_sqrt = np.sqrt(qx**2 + qy**2 + qz**2)

                if q_sqrt > largest_q:
                    continue

                remain = (int(q_sqrt * 1000)) % 100
                if remain >= 90 or remain <= 10:
                    qvalue_10 = int(round(q_sqrt * 10))  # q values multiply by 10

                    if num_qvec[qvalue_10] < maximum_q:                    
                        counter += 1
                        num_qvec[qvalue_10] += 1
    return counter



def sqtgenerator(number_q, cell):
    """
    This function generates the q vector
    """
    total_q = number_q
    minq = 2 * np.pi / cell  # 2pi/L
    largest_q = 2.05 # calculates sqt out to q = 2 Angstroms^-1
    maximum_q = 100 # 100 q vectors per bin of q
    n_max = int(np.ceil(largest_q/minq))

    q = np.zeros(total_q); qtmps = []
    q_vec = np.zeros((total_q, 3)); q_vec_tmps = []
    for nx in range(-n_max, n_max+1, 1):
        for ny in range(-n_max, n_max+1, 1):
            for nz in range(-n_max, n_max+1, 1):

                if nx == 0 and ny == 0 and nz == 0:
                    continue

                qx = nx * minq; qy = ny * minq; qz = nz * minq
                q_sqrt = np.sqrt(qx**2 + qy**2 + qz**2)

                if q_sqrt > largest_q:
                    continue
                remain = (int(q_sqrt * 1000)) % 100  # assure the generated q has a variance of +- 10%
                
                if remain >= 90 or remain <= 10:
                    qvalue_10 = int(round(q_sqrt * 10))  # q values multiply by 10

                    qtmps.append(q_sqrt)
                    q_vec_tmps.append([qx,qy,qz])

    # randomly selects maximum_q number of q_vectors for each bin of q
    qtmps = np.array(qtmps); q_vec_tmps = np.array(q_vec_tmps); count = 0
    for i in range(21):
        arr = np.where((qtmps > 0.1*i-0.05) & (qtmps < 0.1*i+0.05))[0]

        if len(arr) > maximum_q:
            arr = np.array(random.sample(list(arr), k=maximum_q))

        for j in arr:
            q[count] = qtmps[j]
            q_vec[count] = q_vec_tmps[j]
            count += 1
            
    return q_vec, q, total_q



def calc_sqt(dt_max, nframes, scatter_length, qs):
    f = h5py.File('/tmp/sq.hdf5','r'); dset1 = f['master']

    sqt_arr = np.zeros(dt_max)
    cossum = np.zeros(nframes); sinsum = np.zeros(nframes)
    for i in range(nframes):
        pos = dset1[i]

        qr = (qs * pos).sum(axis=1)
        cossum[i] = (scatter_length * np.cos(qr)).sum()
        sinsum[i] = (scatter_length * np.sin(qr)).sum()

    for i in range(dt_max):
        for j in range(0, nframes - i):
            k = i + j
            sqt_arr[i] += (2.0 * (cossum[j] * cossum[k] + sinsum[j] * sinsum[k]))
    # averaging for the statistic enhancement
    for i in range(dt_max):
        sqt_arr[i] /= (2.0 * (nframes - i))
    
    f.close()

    return sqt_arr



def sqt_analysis(dt, dt_max, cell, nframes, scatter_length):
# Calculate the dynamic structure factor
# inputs: dt = timestep; dt_max = simulation time; cell = simulation cell length; nframes = total number of frames; scatter_length = scatter factor of each atom

    # Get q vectors (q = 2pi/L*(nx, ny, nz), where ni is an integrer)
    number_q = presqt(cell)
    q_vec, q, n_max = sqtgenerator(number_q, cell)

    # determined by largest_q, i.e., (5.05 + 0.05)*10
    q_points = 21; dt_max = int(dt_max/dt)
    sqttmps = np.zeros((q_points, dt_max))
    sqtvals = np.zeros((q_points, dt_max))  # bins from 0.0 to 2.0 with increament 0.1

    # calculate sqt at each q_vector
    pool = mp.Pool(processes=nt)
    fn_q_analysis = functools.partial(calc_sqt, dt_max, nframes, scatter_length)
    sqt = pool.map(fn_q_analysis, q_vec); sqt = np.array(sqt)
    pool.close()
    pool.join()

    # bin sqt based on the value of q
    for i in range(q_points):
        for qindex in range(n_max):
            if int(round(q[qindex] * 10)) == i:
                for time in range(dt_max):
                    sqttmps[i][time] += sqt[qindex][time]

    # normallization
    for i in range(len(sqtvals)):
        if sqttmps[i][0] != 0:
            sqtvals[i] = sqttmps[i] / sqttmps[i][0]

    sqtPrint('sqt', sqtvals, dt)



def load_TRR():
# loads in the trajectory and saves the necessary data to a temporary h5py file

    global t_min, t_max, step

    uta = mda.Universe(top_file, trj_file)
    poly = uta.select_atoms(poly_names)

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

    # predefined scatter factor for individual atom (atom mass, scatter length)
    # obtained from https://www.ncnr.nist.gov/resources/n-lengths/list.html
    scatter_factor = dict([(1, -3.7406), (12, 6.6511), # H and C
                           (14, 9.3700), (16, 5.8030), # N and O
                           (19, 5.6540), (7, -1.9000), # F and Li
                           (32, 2.8040)])              # S
    scatter_length = np.array(list(map(scatter_factor.get, poly.masses.round().astype(int))))
    cell = uta.dimensions[0]
    if len(scatter_length[scatter_length==None]) != 0:
        print('Missing one or more atomic scatter factors! Exiting Code')
        exit()

    master = np.zeros((len(frame_ids),len(poly),3))
    for i, frame in enumerate(frame_ids):
        ts = uta.trajectory[frame]

        if frame%1000 == 0:
            print("Frame " + str(frame))
        
        master[i] = poly.positions

    with h5py.File('sq.hdf5','w') as f:
        dset1 = f.create_dataset("dt", data = [dt])
        dset2 = f.create_dataset("frames", data = frame_ids)
        dset3 = f.create_dataset("master", data = master)
        dset4 = f.create_dataset("cell", data = [cell])
        dset5 = f.create_dataset("scatter_length", data = scatter_length)




def main(trj_file, top_file, poly_names, t_min, t_max, step, nt):

    # Makes sure the temporary h5py file does not exist

    # Load in the trajectory file
    if not os.path.exists('sq.hdf5'):
        load_TRR()
        exit()
    
    if not os.path.exists('/tmp/sq.hdf5'):
        file_size = os.path.getsize('sq.hdf5') / (1024**3)
        print(file_size)
        if file_size > 125.0:
            print("File too large: " + file_size + "GB")
            exit()
        else:
            os.system('cp sq.hdf5 /tmp/')
    
    with h5py.File('/tmp/sq.hdf5','r') as f:
        dset1 = f['frames']; frame_ids = dset1[:]
        dset2 = f['dt']; dt = dset2[0]
        dset3 = f['cell']; cell = dset3[0]
        dset4 = f['scatter_length']; scatter_length = dset4[:]
    nframes = len(frame_ids)

    # If static structure factor does not exist, calculate it. Exits after execution to purge memory
    if not os.path.exists('sq.xvg'):
        sq_analysis(cell, nframes, scatter_length)
        exit()
    
    # If dynamic strucutre factor does not exist, calculate it. Exits after execution to purge memory
    if not os.path.exists('sqt.xvg'):
        if t_calc == -1:
            if step == -1:
                dt_max = -(frame_ids[-1] - frame_ids[0])*dt/step # Time to calculate st out to in ps (10000 is 10 ns)
            else:
                dt_max = (frame_ids[-1] - frame_ids[0])*dt/step # Time to calculate st out to in ps (10000 is 10 ns)
        else:
            dt_max = t_calc

        sqt_analysis(dt, dt_max, cell, nframes, scatter_length)

    os.remove('sq.hdf5'); os.remove('/tmp/sq.hdf5')
    
if __name__ == "__main__":
    main(trj_file, top_file, poly_names, t_min, t_max, step, nt)
