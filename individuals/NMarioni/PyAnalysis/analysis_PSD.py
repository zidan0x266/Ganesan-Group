# This script determines the pore size distribution based on the van der Waals volume of the system
# This code was specifically desgined to find the distribution of water-rich pores within a hydrated polymer system
# The output includes the Cumulative Pore Size Distribution, Pore Size Distribution, and Free Volume Fraction
# This code was written based on the methods used for PoreBlazer: https://github.com/SarkisovGitHub/PoreBlazer

import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
import MDAnalysis.lib.distances as distances
import numpy as np
import secrets

import multiprocessing as mp
import functools
import h5py
import os
from sys import argv

script, trj_file, top_file, system_name, probe_name, probe_radius, t_min, t_max, N_frames, nt = argv
probe_radius = float(probe_radius); t_min = float(t_min); t_max = float(t_max); N_frames = int(N_frames); nt = int(nt)
# trj_file = .trr/.xtc file, top_file = .tpr file
# system_name = names of atoms that make up the polymer matrix in a form acceptable by MDAnalysis.select_atoms(), e.g., "moltype MOL"
# probe_name = names of atoms that make up the "probes" within the polymer matrix in a form acceptable by MDAnalysis.select_atoms(), e.g., water molecules + salt as "name OW or name LI or name CL"
# probe_radius = radius in Angstroms of the probe for obtaining the probe-accessible PSD and FFV, e.g., water molecules are 1.4 Angstroms
# t_min = start time for analysis in ps (-1 assumes the start time of the first frame)
# t_max = end time for analysis in ps (-1 assumes the end time of the last frame)
# N_frames = number of frames to analyze (-1 assumes N_frames = nt) -> for efficiency, N_frames should be a multiple of the number of threads
#            N_frames = 1 defaults to the frame at t_max
# nt = number of threads

# Example of Use: python3 analysis_PSD.py md.xtc md.tpr "moltype MOL" "name OW or name LI or name CL" 1.4 90000 100000 80 80
#                 Polymer chains are molecules with name MOL -> this defines the "system"
#                 Water molecules are defined by the water oxygen atoms, and atoms as is -> this defines the number of molecules within the free volume
#                 probe_radius = 1.4 is approximately the size of a water molecule
#                 The run time of this code is highly dependent on the number of frames -> N_frames shouldn't be too many multiples of nt
#                     Note: each frame takes about the same amount of time to process, so it is efficient for the number of frames examined to be a multiple of nt

# Van der waals radii from Bondi (1964)
# https://www.knowledgedoor.com/2/elements_handbook/bondi_van_der_waals_radius.html
Size_arr = np.array([('C',1.7), ('O',1.52), ('H', 1.2), ('N', 1.55), ('S', 1.8), ('F', 1.47), ('P', 1.80)])



def iDist(frame):
    if frame%nt == 0:
        print("Frame " + str(frame))
    
    rng = np.random.default_rng(secrets.randbits(128))                                                  # Create random seed for random number generator per NumPy documentation

    with h5py.File('/tmp/PSD.hdf5','r') as f:
        dset1 = f['system']; sys = dset1[frame]                                                         # Position of all system atoms
        dset2 = f['sys_radii']; sys_radii = dset2[:]                                                    # Radius of all system atoms, where distance between sphere (see below) and polymer atom minus the radius is the distance to the van der waals surface of the atom
        dset3 = f['cells']; cell = dset3[frame]                                                         # Size of the cell
        dset4 = f['probe']; N_probe = len(dset4[frame])                                                 # Number of "probe" atoms within the free volume
    




    # This part of the calculation determines the maximum size of random spheres without overlapping system atoms, where the total volume of all spheres larger than probe_radius defines the probe-accessible free volume of the system
    # This code will generate random points and grow these points into free volume spheres 
    # This code may need to be adapted into a loop for systems with a very large free volume

    N_spheres = 100 * N_probe                                                                           # Create 100 * (# of "probe" atoms) spheres to sample the free volume space
    sphere_arr = np.zeros((N_spheres,3))                                                                # Array containing coordinates of the center of each sphere
    radii_arr = np.zeros((N_spheres))                                                                   # Array containing the maximum radius of each sphere without overlapping system atoms

    ## Fill sphere_arr with random coordinates
    sphere_arr = np.round(rng.random((N_spheres,3))*cell[:3], decimals=5)
    
    # To reduce the number of calculations and limit memory usage, the distance between sphere centers and polymer atoms is done in steps of 2 Angstroms
    d = 0.0                                                                                             # Maximum distance to calculate between every sphere center and every system atom
    remaining = np.where(radii_arr == 0)[0]                                                             # Index of sphere centers that still need their size determined
    while len(remaining) > 0:
        d += 2

        pair_arr, dist_arr = distances.capped_distance(sphere_arr[remaining], sys, d, box=cell)         # Distance between sphere centers and system atoms
        dist_arr -= sys_radii[pair_arr[:,1]]                                                            # Subtract radius of each system atom from the distance to get the distance to the surface of the atom

        ## Useful print command for troubleshooting memory problems: prints the distance calulcated out to, the number of sphere centers, and the total number of distances to system atoms generated
        #if frame == 0:
        #    print("Create Spheres:", d, len(remaining), len(dist_arr))

        # Fill radii_arr for all sphere centers that contain system atoms within d Angstroms
        index = 0; c_index = pair_arr[0,0]; skip = 0
        for i,c in enumerate(pair_arr[:,0]):
            if c > c_index:
                r_min = np.round(np.min(dist_arr[index:i]), decimals=5)                                 # Minimum distance between sphere center and system atom
                if r_min > 0:
                    radii_arr[remaining[c_index]] = r_min                                               # Sphere does not overlap the system
                else:
                    radii_arr[remaining[c_index]] = -1                                                  # Sphere is within the van der waals surface of the system
                index = i; c_index = c

                if i == len(pair_arr[:,0]) - 1:                                                         # Check to make sure the last sphere center is counted
                    skip = 1
        if skip == 0:
            r_min = np.round(np.min(dist_arr[index:]), decimals=5)
            if r_min > 0:
                radii_arr[remaining[c_index]] = r_min                                                   # Minimum distance between sphere center and system atom
            else:
                radii_arr[remaining[c_index]] = -1                                                      # Sphere does not overlap the system
        del pair_arr; del dist_arr
        
        remaining = np.where(radii_arr == 0)[0]                                                         # Sphere is within the van der waals surface of the system
    del remaining

    # Reduce sphere_arr to spheres of the desirable size, typically minimum radius of probe_radius
    sphere_arr = sphere_arr[np.where(radii_arr > probe_radius)[0]]; radii_arr = radii_arr[radii_arr > probe_radius]; max_radius = np.max(radii_arr)
    # Useful print command for troubleshooting problems: prints the number of sphere_centers within the free volume, the radius of the largest sphere, and the diameter of the largest sphere (pore)
    if frame%nt == 0:
        print("Spheres Created:", len(radii_arr), max_radius, 2 * max_radius)

    ### Code to write coordinates and radius of each random sphere to a .xyz file, which can be visualized in Ovito
    #if frame == 0:
    #    with open('test.xyz', 'w') as anaout:
    #        print(str(len(sphere_arr)), file=anaout)
    #        print('Properties=species:S:1:pos:R:3:Radius:R:1', file=anaout)
    #        for i, sph in enumerate(sphere_arr):
    #            print('X {:10.5f} {:10.5f} {:10.5f} {:10.5f}'.format(sph[0] - cell[0]/2, sph[1] - cell[1]/2, sph[2] - cell[2]/2, radii_arr[i]), file=anaout)
    #    print('XYZ File Printed')





    # This part of the calculation determines the probe-accessible fractional free volume
    # The code will create (number of system + "probe" atoms) random coordinates and determine if they lie within the free volume or not
    # This code is carried out in batches of size N_FFV, which can be adjusted depending on memory usage

    N_FFV = int((len(sys) + N_probe)/2); FFV_track = 0; FFV_total = 0
    while FFV_total < 10 * (len(sys) + N_probe):
        FFV_check = np.zeros((N_FFV), dtype=bool); FFV_total += N_FFV                                                               # FFV_check == False if this coordinate has not been evaluated yet; FFV_total tracks the total number of coordinates generated
        FFV_probe = np.round(rng.random((N_FFV,3))*cell[:3], decimals=5)                                                            # Generate random coordinates within the system volume

        # To reduce the number of calculations and limit memory usage, random coordinates within the van der waals surface of the system are eliminated first
        pair_arr, dist_arr = distances.capped_distance(FFV_probe, sys, np.max(np.array(Size_arr[:,1], dtype=float))+0.5, box=cell)  # Distance between each random coordinate and the system atoms
        dist_arr -= sys_radii[pair_arr[:,1]]                                                                                        # Subtract radius of each system atom to find distance between random coordinates and the van der waals surface of the atoms
        ## Useful print command for troubleshooting memory problems: prints the maximum distance calculated between random coordinates and system atoms, the number of random coordinates, and the number of distances generated
        #if frame == 0:
        #    print("FFV Probes not in System:", np.max(np.array(Size_arr[:,1], dtype=float))+0.5, len(FFV_probe), len(dist_arr))
    
        # Fill FFV_check for all random coordinates located within the van der waals surface of the system
        index = 0; c_index = pair_arr[0,0]; skip = 0
        for i,c in enumerate(pair_arr[:,0]):
            if c > c_index:
                if np.any(dist_arr[index:i] < 0):                                                                                   # If this coordinate is within the van der waals surface of any system atom, mark it as no longer needing to be calculated
                    FFV_check[c_index] = True
                index = i; c_index = c
    
                if i == len(pair_arr[:,0]) - 1:                                                                                     # Check to make sure the last sphere center is counted
                    skip = 1
        if skip == 0:
            if np.any(dist_arr[index:] < 0):                                                                                        # If this coordinate is within the van der waals surface of any system atom, mark it as no longer needing to be calculated
                FFV_check[c_index] = True
        del pair_arr; del dist_arr
    
        # To reduce the number of calculations and limit memory usage, the distance between sphere centers and polymer atoms is done in steps of 2 Angstroms
        d = 0                                                                                                                       # Maximum distance to calculate between every coordinate and every sphere center
        remaining = np.where(FFV_check == False)[0]                                                                                 # Index of coordinates that still need their location determined
        r_OLD = -1                                                                                                                  # Track how many points are remaining at previous value of d
        while d < max_radius:
            d += 2

            # Some number of spheres will never be assigned due to residing in free volume elements smaller than probe_radius
            # Setting d to probe_radius + 0.5 once most of the points have been assigned dramatically speeds up the calculations
            # The value checked against, here 0.05, may need to be adjusted depending on the system
            if r_OLD > 0 and d > 2*probe_radius and (r_OLD - len(remaining))/len(remaining) < 0.1:
                d = max_radius + 0.5
    
            pair_arr, dist_arr = distances.capped_distance(FFV_probe[remaining], sphere_arr, d, box=cell)                          # Distance between each random coordinate and the free volume sphere centers generated in the code above
            dist_arr -= radii_arr[pair_arr[:,1]]                                                                                   # Subtract radius of each sphere to find distance between random coordinate and the surface of the sphere
            ## Useful print command for troubleshooting memory problems: prints the maximum distance calculated between random coordinates and system atoms, the number of random coordinates, and the number of distances generated
            #if frame%nt == 0:
            #    print("FFV Probes:", d, len(remaining), len(dist_arr))
    
            # Fill FFV_check for all random coordinates located within a free volume sphere
            index = 0; c_index = pair_arr[0,0]; skip = 0
            for i,c in enumerate(pair_arr[:,0]):
                if c > c_index:
                    if np.any(dist_arr[index:i] < 0):                                                                               # If the random coordinate lies within a free volume sphere, mark it as no longer needing to be calculated and increase the number of coordinates within the free volume by 1
                        FFV_check[remaining[c_index]] = True
                        FFV_track += 1
                    index = i; c_index = c
    
                    if i == len(pair_arr[:,0]) - 1:                                                                                 # Check to make sure the last coordinate is counted
                        skip = 1
            if skip == 0:
                if np.any(dist_arr[index:] < 0):                                                                                    # If the random coordinate lies within a free volume sphere, mark it as no longer needing to be calculated and increase the number of coordinates within the free volume by 1
                    FFV_check[remaining[c_index]] = True
            del pair_arr; del dist_arr
    
            r_OLD = len(remaining); remaining = np.where(FFV_check == False)[0]
        ## Useful print command to track the probe-accessible free volume every loop. The maximum number of coordinates of N_FFV can be altered as needed to balance speed and accuracy
        #if frame == 0:
        #    print("FFV:", FFV_total, FFV_track / FFV_total)
    # Useful print command to track the final FFV
    if frame%nt == 0:
        print("FFV Final:", FFV_total, FFV_track / FFV_total)
    del remaining; del FFV_check





    # This part of the calculation determines the cumulative probe-accessible pore size distribution, where the distribution is defined as the probability that a random point within the free volume resides within a free volume sphere of diameter d with minimum size probe_radius
    # This code will generate (10 * N_probe) random points within the free volume and determine the largest free volume sphere that contains that point
    # This code is carried out in batches of size N_PSD, which can be adjusted depending on memory usage

    d_arr = np.arange(0, 50.25, 0.25); PSD_arr = np.zeros_like(d_arr, dtype=int)                                                    # d_arr is the histogram of free volume sphere sizes; PSD_arr tracks the number of instances of points contained within free volume spheres of size at least d

    N_PSD = int(N_probe/2)
    while PSD_arr[0] < 10 * N_probe:
        PSD_probes = np.zeros((N_PSD,3))                                                                                            # Array containing coordinates of each random point

        # Fill PSD_probes with random coordinates within probe_radius of at least 1 sphere
        index = 0
        while np.all(PSD_probes[-1] == 0):
            PSD_probes_temp = np.round(rng.random((N_PSD - index,3))*cell[:3], decimals=5)                                          # Generate random coordinates within the system volume
        
            pair_arr, dist_arr = distances.capped_distance(PSD_probes_temp, sphere_arr, probe_radius, box=cell)                     # Find distance between random coordinates and free volume spheres
            PSD_probes_temp = PSD_probes_temp[np.unique(pair_arr[:,0])]                                                             # Eliminate coordinates not within probe_radius of at least 1 sphere
        
            PSD_probes[index:index + len(PSD_probes_temp)] = PSD_probes_temp; index += len(PSD_probes_temp)                         # Update PSD_probes and loop around until PSD_probes is full
        del PSD_probes_temp; del index; del pair_arr; del dist_arr

        pair_arr, dist_arr = distances.capped_distance(PSD_probes, sphere_arr, max_radius+0.5, box=cell)                            # Find distance between random coordinates and free volume sphere centers
        dist_arr -= radii_arr[pair_arr[:,1]]                                                                                        # Subtract radius of each free volume sphere to find distance between random coordinates and the surface of the sphere
        ## Useful print command for troubleshooting memory problems: prints the maximum distance calculated between random coordinates and free volume sphere centers, the number of random coordinates, and the number of distances generated
        #if frame == 0:
        #    print("PSD probes:", max_radius+0.5, N_PSD, len(dist_arr))

        # Fill PSD_arr for all random coordinates located within a free volume sphere
        index = 0; c_index = pair_arr[0,0]; skip = 0
        for i,c in enumerate(pair_arr[:,0]):
            if c > c_index:
                contained = pair_arr[np.where(dist_arr[index:i] < 0)[0],1]                                                          # Index of all free volume spheres containing the points
                if len(contained) > 0:                                                                                              # If at least 1 free volume sphere contains the point, find the largest sphere and add it to PSD_arr in a cumulative manner
                    max_size = 2 * np.max(radii_arr[contained])
                    PSD_arr[np.where(d_arr <= max_size)[0]] += 1
                index = i; c_index = c

                if i == len(pair_arr[:,0]) - 1:                                                                                     # Check to make sure the last coordinate is counted
                    skip = 1
        if skip == 0:
            contained = pair_arr[np.where(dist_arr[index:] < 0)[0],1]                                                               # Index of all free volume spheres containing the points
            if len(contained) > 0:                                                                                                  # If at least 1 free volume sphere contains the point, find the largest sphere and add it to PSD_arr in a cumulative manner
                max_size = 2 * np.max(radii_arr[contained])
                PSD_arr[np.where(d_arr <= max_size)[0]] += 1
        del pair_arr; del dist_arr

        # Useful print command to track the probe-accessible PSD every loop. The maximum number of coordinates of N_PSD can be altered as needed to balance speed and accuracy
        if frame == 0:
            print("PSD:", PSD_arr[0])
        #    print_string=''
        #    for i in PSD_arr:
        #        if i == 0:
        #            continue
        #        print_string += str(np.round(i / PSD_arr[0], decimals=5)) + ' '
        #    print(print_string)
    # Useful print command to track the final PSD
    if frame%nt == 0:
        print("PSD Final:", PSD_arr[0])
    #    print_string=''
    #    for i in PSD_arr:
    #        if i == 0:
    #            continue
    #        print_string += str(np.round(i / PSD_arr[0], decimals=5)) + ' '
    #    print(print_string)
    del contained




    # Return the necessary information to complete the calculations: FFV_track / FFV_total gives the pore-accessible free volume, PSD_arr / PSD_arr[0] gives the pore-accessible PSD
    return np.hstack((np.array([FFV_track, FFV_total], dtype=int),PSD_arr), dtype=int)
    


def load_TRR():
# loads in the trajectory and saves the necessary data to a temporary h5py file

    global t_min, t_max, N_frames

    uta = mda.Universe(top_file, trj_file)  # Load in the trajectory and topology
    probe = uta.select_atoms(probe_name)    # Define the "probe" atoms
    system = uta.select_atoms(system_name)  # Define the system atoms
    
    # Create an array that tracks the radius of each system atom based on Size_array
    sys_names = system.names; sys_radii = []
    for name in sys_names:
        name = str(name)[0]
        if name in Size_arr[:,0]:
            sys_radii.append(float(Size_arr[np.where(Size_arr[:,0] == name)[0][0],1]))
        else:
            print("Missing Atom Name and Size in Size_arr (See Atom Name Below)")
            print(name)
            exit()
    sys_radii = np.array(sys_radii)

    # Define the system times/frames to be calculated over
    if t_min == -1:
        t_min = uta.trajectory[0].time
    if t_max == -1:
        t_max = uta.trajectory[-1].time
    if N_frames < 1:
        N_frames = nt
    dt = np.round((uta.trajectory[1].time - uta.trajectory[0].time),3)
    if N_frames == 1:
        frame_ids = np.array([int((t_max - uta.trajectory[0].time)/dt)], dtype=int)
    else:
        frame_ids = np.linspace(int((t_min - uta.trajectory[0].time)/dt), int((t_max - uta.trajectory[0].time)/dt), N_frames, dtype=int)
        print("Timestep: ~" + str(dt*(frame_ids[1] - frame_ids[0])) + " ps")
    print("Number of Frames: " + str(len(frame_ids)))

    # Load in the necessary data: "probe" atom position, system atom positions, cell dimensions
    r_probe = []; r_system = []; cells = []
    for frame in frame_ids:

        if frame%5000 == 0:
            print("Frame " + str(frame))

        ts = uta.trajectory[frame]
        cell = ts.dimensions
        cells.append(cell)

        r_probe.append(probe.positions)
        r_system.append(system.positions)

    # Save necessary infomration to a .hdf5 file for later use in the calculation
    with h5py.File('/tmp/PSD.hdf5','w') as f:
        dset1 = f.create_dataset("probe", data=r_probe)
        dset2 = f.create_dataset("system", data=r_system)
        dset3 = f.create_dataset("cells", data = cells)
        dset4 = f.create_dataset("frames", data = frame_ids)
        dset5 = f.create_dataset("sys_radii", data = sys_radii)



def main(trj_file, top_file, system_name, probe_name, probe_radius, t_min, t_max, N_frames, nt):

    # Load in the trajectory file and exit the code to purge the memory before multi-threading
    if not os.path.exists('/tmp/PSD.hdf5'):
        load_TRR()
        exit()

    with h5py.File('/tmp/PSD.hdf5','r') as f:
        dset1 = f['frames']; frame_ids = dset1[:]
    frame_ids = np.arange(0,len(frame_ids),1)

    print("FFV/PSD Analysis")
    # Perform the analysis using multi-threading
    pool = mp.Pool(processes=nt)
    func = functools.partial(iDist)
    radii_arr = pool.map(func, list(frame_ids))
    pool.close()
    pool.join()
    radii_arr = np.array(radii_arr)

    # Return the average and standard deviation (over the frames processed) of the probe-accessible fraction free volume
    FFV = radii_arr[:,:2]; FFV = FFV[:,0] / FFV[:,1]; FFV = np.array([np.mean(FFV), np.std(FFV)])
    with open('FFV.xvg', 'w') as anaout:
        print("# FFV Std", file=anaout)
        print('0.0 {:10.5f} {:10.5f}'.format(FFV[0], FFV[1]), file=anaout)

    # Return the average and standard deviation (over the frames processed) of the probe-accessible pore size ditribution
    d_arr = np.arange(0, 50.25, 0.25)
    PSD_all = radii_arr[:,2:]; PSD_all = np.divide(PSD_all.T, PSD_all[:,0], dtype=float).T
    PSD_Cumulative = np.array([np.mean(PSD_all, axis=0), np.std(PSD_all, axis = 0)])
    PSD = np.array([np.mean((PSD_all[:,:len(d_arr)-2] - PSD_all[:,2:])/(d_arr[2:] - d_arr[:len(d_arr)-2]), axis=0), np.std((PSD_all[:,:len(d_arr)-2] - PSD_all[:,2:])/(d_arr[2:] - d_arr[:len(d_arr)-2]), axis=0)])

    with open('Cumulative_PSD.xvg', 'w') as anaout:
        print("# Cumulative_PSD Std", file=anaout)
        for i in range(len(d_arr)):
            print(' {:10.5f} {:10.5f} {:10.5f}'.format(np.round(d_arr[i], decimals=3), PSD_Cumulative[0,i], PSD_Cumulative[1,i]), file=anaout)

    with open('PSD.xvg', 'w') as anaout:
        print("# PSD Std", file=anaout)
        for i in range(len(d_arr)-1):
            if i == 0:
                print(' {:10.5f} {:10.5f} {:10.5f}'.format(d_arr[i], 0.0, 0.0), file=anaout)
            else:
                print(' {:10.5f} {:10.5f} {:10.5f}'.format(np.round(d_arr[i], decimals=3), PSD[0,i-1], PSD[1,i-1]), file=anaout)

    # Deletes the temporary h5py file
    os.remove('/tmp/PSD.hdf5')

if __name__ == "__main__":
    main(trj_file, top_file, system_name, probe_name, probe_radius, t_min, t_max, N_frames, nt)
