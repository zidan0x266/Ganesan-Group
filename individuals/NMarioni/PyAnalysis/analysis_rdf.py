#This script calculates the Radial Distribution Function (RDF) and Cumulative/Coordination Distribution Functions (CDF)

import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
import MDAnalysis.lib.distances as distances
import numpy as np
import h5py
import os
import multiprocessing as mp
import functools

from sys import argv

script, trj_file, top_file, a1_name, a2_name, rmax, cdf_check, t_min, t_max, step, nt = argv
rmax = float(rmax)*10; rlen = int(rmax/0.02); cdf_check = int(cdf_check); t_min = float(t_min); t_max = float(t_max); step = int(step); nt = int(nt)
# trj_file = .trr/.xtc file, top_file = .tpr file
# a1_name = name of reference atom to be analyzed
# a2_name = name of selection atom to be analyzed
# rmax = max distanced to be analyzed in nm
# cdf_check = 0 for no cdf, 1 for coordination of a2 around a1, 2 for coordination of a1 around a2, 3 for both
# (cdf = coordination/cumulative disrubtion function)
# t_min = start time for analysis in ps (-1 assumes the start time of the first frame)
# t_max = end time for analysis in ps (-1 assumes the end time of the last frame)
# step = frame step-size (-1 assumes a step-size of 1)
# nt = number of threads

# Example: python3 ${Path}/analysis_rdf.py md.trr md.tpr NA CL 2.00 1 180000 200000 1 128



def rdfPrint(filename, bins, rdf):  # print displacement
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# r g(r)', file=anaout)
        for i in range(0, len(bins)):
            print('{:10.3f} {:10.5f}'.format(bins[i], rdf[i]), file=anaout)



def cdfPrint(filename, bins, rdf):  # print displacement
    with open('{}.xvg'.format(filename), 'w') as anaout:
         print('# r n(r)', file=anaout)
         for i in range(0, len(bins)):
            print('{:10.3f} {:10.5f}'.format(bins[i], rdf[i]), file=anaout)



def rdf_analysis(nframe):
# calculate the radial distribution functions
#
# Inputs: nframe = number of frames being analyzed

    rdfrange = (0.0, rmax)
    edges = np.histogram(rdfrange, rlen)[1]
    bins = 0.5 * (edges[:-1] + edges[1:])
    dr = edges[1]
    vol = np.power(edges[1:], 3) - np.power(edges[:-1], 3)
    vol *= 4/3.0 * np.pi

    #with h5py.File('/tmp/r_RDF_'+a1_name+'_'+a2_name+'.hdf5','r') as f:
    #    dset1 = f['r1']; r1 = dset1[:]
    #    dset2 = f['r2']; r2 = dset2[:]
    #    dset3 = f['cells']; cells = dset3[:]

    pool = mp.Pool(processes=nt)
    func = functools.partial(rdf_calc, vol)
    rdf = pool.map(func, list(range(nframe)))
    pool.close()
    pool.join()

    rdf = np.average(rdf, axis = 0)
    rdf_name = 'rdf_{}_{}'.format(a1_name, a2_name)
    rdfPrint(rdf_name, np.divide(bins, 10), rdf)

    # Calculate the CDF
    if cdf_check != 0:
        with h5py.File('/tmp/r_RDF_'+a1_name+'_'+a2_name+'.hdf5','r') as f:
            dset1 = f['r1']; r1 = dset1[:]
            dset2 = f['r2']; r2 = dset2[:]
            dset3 = f['cells']; cells = dset3[:]

            # Retrieve the average number density of a1 and a2
            a1_rho = []; a2_rho = []
            for cell in cells:
                tot_vol = cell[0] * cell[1] * cell[2]
                a1_rho.append(len(r1[0,:,0]) / tot_vol)
                a2_rho.append(len(r2[0,:,0]) / tot_vol)
            a1_rho = np.mean(a1_rho); a2_rho = np.mean(a2_rho)
        
        cdf_tmp = rdf
        cdf_tmp = np.multiply(cdf_tmp, np.power(bins, 2))
        cdf_tmp = np.multiply((cdf_tmp[1:] + cdf_tmp[:-1])/2, bins[1:] - bins[:-1])
        cdf_tmp = np.insert(cdf_tmp, 0, 0.0)
        cdf_tmp = np.cumsum(cdf_tmp)
        cdf_tmp *= 4*np.pi
        if(cdf_check == 1 or cdf_check == 3):
            cdf1 = (np.multiply(cdf_tmp, a2_rho))
            cdf_name = 'cdf_{}_{}'.format(a1_name, a2_name)
            cdfPrint(cdf_name, np.divide(bins, 10), cdf1)
        if(cdf_check == 2 or cdf_check == 3):
            cdf2 = (np.multiply(cdf_tmp, a1_rho))
            cdf_name = 'cdf_{}_{}'.format(a2_name, a1_name)
            cdfPrint(cdf_name, np.divide(bins, 10), cdf2)



def rdf_calc(vol, frame):
# calculate the radial distribution function
#
# Inputs: vol = array of volumes for each bin; frame = frame being analyzed

    if frame%500 == 0:
            print("Frame "+ str(frame))
    
    with h5py.File('/tmp/r_RDF_'+a1_name+'_'+a2_name+'.hdf5','r') as RDF_file:
        r1 = RDF_file['r1'][frame,:,:]; r2 = RDF_file['r2'][frame,:,:]; cell = RDF_file['cells'][frame]

    tot_vol = cell[0] * cell[1] * cell[2]
    rho = len(r1) * len(r2) / tot_vol

    tmp = []
    for a1_i in r1:
        a2_tmp = r2
        if a1_name == a2_name:
            a2_tmp = np.delete(a2_tmp, np.where((a2_tmp == a1_i).all(axis=1))[0][0], 0)
        distpair = distances.capped_distance(a1_i, a2_tmp, rmax, box=cell)[1]
        count = np.histogram(distpair, rlen, (0, rmax))[0]
        tmp.append(np.divide(count,np.multiply(vol,rho)))

    return np.sum(tmp, axis = 0)



def load_TRR():
# Load in the trajectory file and write necessary data to a h5py file
# For very large files, it is recommended to dump only the atoms you are interested in to a separate .trr/.xtc file for analysis

    global t_min, t_max, step

    uta = mda.Universe(top_file, trj_file, tpr_resid_from_one=True)
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
    
    r1_temp = []; r2_temp = []; cells = []
    for frame in frame_ids:
        ts = uta.trajectory[frame]
        cell = ts.dimensions
        cells.append(cell)
    
        if ts.time%5000 == 0:
            print("Time "+ str(ts.time))

        r1_temp.append(a1.positions)
        
        r2_temp.append(a2.positions)

    r1 = np.array(r1_temp)
    r2 = np.array(r2_temp)

    with h5py.File('/tmp/r_RDF_'+a1_name+'_'+a2_name+'.hdf5','w') as f:
        dset1 = f.create_dataset("r1", data=r1)
        dset2 = f.create_dataset("r2", data=r2)
        dset3 = f.create_dataset("dt", data=[dt])
        dset4 = f.create_dataset("frames", data=frame_ids)
        dset5 = f.create_dataset("cells", data = cells)



def main(trj_file, top_file, a1_name, a2_name, rmax, cdf_check, t_min, t_max, step, nt):

    # If there is not a position h5py file, then create one and end the program
    # This is done to avoid memory problems during multiprocessing
    if not os.path.exists('/tmp/r_RDF_'+a1_name+'_'+a2_name+'.hdf5'):
        load_TRR()
        #exit()
    with h5py.File('/tmp/r_RDF_'+a1_name+'_'+a2_name+'.hdf5','r') as f:
        dset2 = f['frames']; frame_ids = dset2[:]
    
    print("RDF Analysis")
    nframe = len(frame_ids)
    rdf_analysis(nframe)

    os.remove('/tmp/r_RDF_'+a1_name+'_'+a2_name+'.hdf5')

if __name__ == "__main__":
    main(trj_file, top_file, a1_name, a2_name, rmax, cdf_check, t_min, t_max, step, nt)
