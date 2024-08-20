"""
Copyright (C) 2018-2021 Zidan Zhang <zhangzidan@gmail.com>

This file is distributed under free software licence:
you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# Original Code can be found here: https://github.com/kdfong/transport-coefficients-MSD
# Adjustments were made to load data similarly to my other codes. The base FFT methods have not been changed.
# This script calculates "MSDs" for retrieving Onsager Coefficients

import multiprocessing as mp
import functools

import numpy as np
import h5py
import MDAnalysis as mda
import MDAnalysis.lib.distances as distances
import os

from sys import argv
script, trj_file, top_file, cation_name, anion_name, polyanion_name, t_min, t_max, step, nt, CoM_check = argv
t_min = float(t_min); t_max = float(t_max); step = int(step); nt = int(nt)
# trj_file = .trr/.xtc file, top_file = .tpr file
# cation_name = name of cation
# anion_name = name of anion
# polyanion_name = name of polyanion
# t_min = start time for analysis in ps (-1 assumes the start time of the first frame)
# t_max = end time for analysis in ps (-1 assumes the end time of the last frame)
# step = frame step-size (-1 assumes a step-size of 1)
# nt = number of threads
# NOTE: Make sure to have unwrapped coordinates (gmx trjconv with -pbc nojump)
# CoM_check = 0 if the trj_file contains ALL atoms, 1 else
#   #NOTE: If the trj_file does NOT contain ALL atoms, then the system center of mass will not be correct. First run Sys_CoM.py and then send CoM_check = 1

# Example: python3 ${Path}/analysis_transport_3ions.py unwrap.trr md.tpr -1 -1 -1 128 0

def outPrint(filename, dt, msd):  # print msd
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# time msd', file=anaout)
        for i in range(0, len(dt)):
            print(' {: 1.5e} {: 1.5e}'.format(dt[i], msd[i]), file=anaout)


def create_mda(top, trj):
    """
    Creates MDAnalysis universe with trajectory data.
    :param top: string, path to Gromacs topology file with atom coordinates and topology
    :param trj: string or list[string], path(s) to Gromacs xtc files with trajectory data
    :return uta: MDAnalysis universe
    """
    uta = mda.Universe(top, trj)
    return uta


def define_atom_types(uta, cation_type = 'name '+cation_name, anion_type = 'name '+anion_name, polyanion_type = 'name '+polyanion_name):
    """ 
    Sorts atoms in the MDAnalysis universe based on type (cation and anion).
    Selections must be a single atom (rather than a molecule center of mass).
    :param uta: MDAnalysis universe
    :param cation_type: string, ion name corresponding to cations in the Gromacs topology files
    :param anion_type: string, ion name corresponding to anions in the Gromacs topology files
    :param polyanion_type: string, ion name corresponding to polyanions in the Gromacs topology files
    :return cations, anions: MDAnalysis AtomGroups corresponding to cations and anions
    """
    cations = uta.select_atoms(cation_type)
    anions = uta.select_atoms(anion_type)
    polyanions = uta.select_atoms(polyanion_type)
    return cations, anions, polyanions


# Algorithms in this section are adapted from DOI: 10.1051/sfn/201112010 and
# https://stackoverflow.com/questions/34222272/computing-mean-square-displacement-using-python-and-fft

def autocorrFFT(x):
    """
    Calculates the autocorrelation function using the fast Fourier transform.
    
    :param x: array[float], function on which to compute autocorrelation function
    :return: acf: array[float], autocorrelation function
    """
    N = len(x)
    F = np.fft.fft(x, n = 2 * N)  
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res = (res[:N]).real   
    n = N * np.ones(N) - np.arange(0, N) 
    acf = res / n
    return acf


def msd_fft(r):
    """
    Computes mean square displacement using the fast Fourier transform.
    
    :param r: array[float], atom positions over time
    :return: msd: array[float], mean-squared displacement over time
    """
    N = len(r)
    D = np.square(r).sum(axis = 1) 
    D = np.append(D, 0) 
    S2 = sum([autocorrFFT(r[:, i]) for i in range(r.shape[1])])
    Q = 2 * D.sum()
    S1 = np.zeros(N)
    for m in range(N):
        Q = Q - D[m - 1] - D[N - m]
        S1[m] = Q / (N - m)
    msd = S1 - 2 * S2
    return msd

def cross_corr(x, y):
    """
    Calculates cross-correlation function of x and y using the 
    fast Fourier transform.
    :param x: array[float], data set 1
    :param y: array[float], data set 2
    :return: cf: array[float], cross-correlation function
    """
    N = len(x)
    F1 = np.fft.fft(x, n = 2**(N * 2 - 1).bit_length())
    F2 = np.fft.fft(y, n = 2**(N * 2 - 1).bit_length())
    PSD = F1 * F2.conjugate()
    res = np.fft.ifft(PSD)
    res = (res[:N]).real   
    n = N * np.ones(N) - np.arange(0, N)
    cf = res / n
    return cf

def msd_fft_cross(r, k):
    """
    Calculates "MSD" (cross-correlations) using the fast Fourier transform.
    :param r: array[float], positions of atom type 1 over time
    :param k: array[float], positions of atom type 2 over time
    :return: msd: array[float], "MSD" over time
    """
    N = len(r)
    D = np.multiply(r,k).sum(axis=1) 
    D = np.append(D,0) 
    S2 = sum([cross_corr(r[:, i], k[:,i]) for i in range(r.shape[1])])
    S3 = sum([cross_corr(k[:, i], r[:,i]) for i in range(k.shape[1])])
    Q = 2 * D.sum()
    S1 = np.zeros(N)
    for m in range(N):
        Q = Q - D[m - 1] - D[N - m]
        S1[m] = Q / (N - m)
    msd = S1 - S2 - S3
    return msd


def create_position_arrays(uta, cations, anions, polyanions, frame_ids):
    """
    Creates an array containing the positions of all cations and anions over time.
    :param uta: MDAnalysis universe
    :param anions: MDAnalysis AtomGroup containing all anions (assumes anions are single atoms)
    :param cations: MDAnalysis AtomGroup containing all cations (assumes cations are single atoms)
    :param frame_ids: array[int], frames at which position data was collected in the simulation
    :return anion_positions, cation_positions: array[float,float,float], array with all
    cation/anion positions. Indices correspond to time, ion index, and spatial dimension
    (x,y,z), respectively
    """
    if CoM_check == '1':
        print('CoM.hdf5 is used')
        with h5py.File('CoM.hdf5','r') as f:
            times = f['times'][:]
            CoM_ar = f['CoM'][:,:]

    cation_positions = np.zeros((len(frame_ids), len(cations), 3))
    anion_positions = np.zeros((len(frame_ids), len(anions), 3))
    polyanion_positions = np.zeros((len(frame_ids), len(polyanions), 3))

    for i, frame in enumerate(frame_ids):
        if frame%5000 == 0:
            print("Frame " + str(frame))

        ts = uta.trajectory[frame]
        cell = ts.dimensions

        if CoM_check == '0':
            CoM = uta.atoms.center_of_mass()
        else:
            CoM = CoM_ar[np.where(times == ts.time)]

        cation_positions[i, :, :] = (cations.positions - CoM)/10.0
        anion_positions[i, :, :] = (anions.positions - CoM)/10.0
        polyanion_positions[i, :, :] = (polyanions.positions - CoM)/10.0
    
    return cation_positions, anion_positions, polyanion_positions

def calc_Lii_self(atom_positions, times):
    """ 
    Calculates the "MSD" for the self component for a diagonal transport coefficient (L^{ii}).
    :param atom_positions: array[float,float,float], position of each atom over time.
    Indices correspond to time, ion index, and spatial dimension (x,y,z), respectively.
    :param times: array[float], times at which position data was collected in the simulation
    :return msd: array[float], "MSD" corresponding to the L^{ii}_{self} transport 
    coefficient at each time
    """
    Lii_self = np.zeros(len(times))
    n_atoms = np.shape(atom_positions)[1]
    for atom_num in (range(n_atoms)):
        r = atom_positions[:, atom_num, :]
        msd_temp = msd_fft(np.array(r))
        Lii_self += msd_temp
    msd = np.array(Lii_self)
    return msd


def calc_Lii(atom_positions, times):
    """ 
    Calculates the "MSD" for the diagonal transport coefficient L^{ii}. 
    :param atom_positions: array[float,float,float], position of each atom over time.
    Indices correspond to time, ion index, and spatial dimension (x,y,z), respectively.
    :param times: array[float], times at which position data was collected in the simulation
    :return msd: array[float], "MSD" corresponding to the L^{ii} transport 
    coefficient at each time
    """
    r_sum = np.sum(atom_positions, axis = 1)
    msd = msd_fft(r_sum)
    return np.array(msd)


def compute_all_Lij(cation_positions, anion_positions, polyanion_positions, times):
    """
    Computes the "MSDs" for all transport coefficients.
    :param cation_positions, anion_positions: array[float,float,float], position of each 
    atom (anion or cation, respectively) over time. Indices correspond to time, ion index,
    and spatial dimension (x,y,z), respectively.
    :param times: array[float], times at which position data was collected in the simulation
    :param volume: float, volume of simulation box
    :return msds_all: list[array[float]], the "MSDs" corresponding to each transport coefficient,
    L^{++}, L^{++}_{self}, L^{--}, L^{--}_{self}, L^{+-}
    """
    print("MSD Self")
    #msd_self_cation = calc_Lii_self(cation_positions, times)
    #msd_self_anion =  calc_Lii_self(anion_positions, times)
    msd_self_polyanion =  calc_Lii_self(polyanion_positions, times)

    print("MSD Distinct")
    #msd_cation = calc_Lii(cation_positions, times)
    #msd_anion = calc_Lii(anion_positions, times)
    msd_polyanion = calc_Lii(polyanion_positions, times)

    print("MSD Cross")
    #msd_cross_catan = calc_Lij(cation_positions, anion_positions, times)
    msd_cross_catpolyan = calc_Lij(cation_positions, polyanion_positions, times)
    msd_cross_anpolyan = calc_Lij(anion_positions, polyanion_positions, times)

    msd_self_cation=msd_self_polyanion*0; msd_self_anion=msd_self_polyanion*0
    msd_cation=msd_polyanion*0; msd_anion=msd_polyanion*0
    msd_cross_catan=msd_cross_catpolyan*0

    msds_all = [msd_cation, msd_self_cation, msd_anion, msd_self_anion, msd_polyanion, msd_self_polyanion, msd_cross_catan, msd_cross_catpolyan, msd_cross_anpolyan]
    return msds_all


def calc_Lij(a1_positions, a2_positions, times):
    """
    Calculates the "MSD" for the off-diagonal transport coefficient L^{ij}, i \neq j.
    :param cation_positions, anion_positions: array[float,float,float], position of each 
    atom (anion or cation, respectively) over time. Indices correspond to time, ion index,
    and spatial dimension (x,y,z), respectively.
    :param times: array[float], times at which position data was collected in the simulation
    :return msd: array[float], "MSD" corresponding to the L^{ij} transport coefficient at 
    each time.
    """
    r_a1 = np.sum(a1_positions, axis = 1)
    r_a2 = np.sum(a2_positions, axis = 1)
    msd = msd_fft_cross(np.array(r_a1), np.array(r_a2))
    return np.array(msd)


def main(trj_file, top_file, t_min, t_max, step, nt):

    uta = create_mda(top_file, trj_file)
    cations, anions, polyanions = define_atom_types(uta)

    if t_min == -1:
        t_min = uta.trajectory[0].time
    if t_max == -1:
        t_max = uta.trajectory[-1].time
    if step < 1:
        step = 1
    dt = np.round((uta.trajectory[1].time - uta.trajectory[0].time),2)
    print("Timestep " + str(dt))
    frame_ids = np.arange(int((t_min - uta.trajectory[0].time)/dt), int((t_max - uta.trajectory[0].time)/dt + 1), step)
    dt = step*dt
    times = frame_ids*dt

    print('Retrieve Ion Data')
    cation_positions, anion_positions, polyanion_positions = create_position_arrays(uta, cations, anions, polyanions, frame_ids)

    print('Compute Lijs')
    msds_all = compute_all_Lij(cation_positions, anion_positions, polyanion_positions, times)

    #z_a1 = 1; z_a2 = -1; z_a3 = -1
    #v_a1 = 1; v_a2 = 1; v_a3 = -1
    #N_a1 = len(cations); N_a2 = len(anions); N_a3 = len(polyanions)
    #N_a1a2 = ((-z_a1*N_a1 * z_a2*N_a2) / (z_a1*z_a1*N_a1 + z_a2*z_a2*N_a2)) * ((v_a1 + v_a2) / (v_a1 * v_a2))
    #N_a1a3 = ((-z_a1*N_a1 * z_a3*N_a3) / (z_a1*z_a1*N_a1 + z_a3*z_a3*N_a3)) * ((v_a1 + v_a3) / (v_a1 * v_a3))
    #N_a2a3 = ((-z_a2*N_a2 * z_a3*N_a3) / (z_a2*z_a2*N_a2 + z_a3*z_a3*N_a3)) * ((v_a2 + v_a3) / (v_a2 * v_a3))

    N_a1 = 1; N_a2 = 1; N_a3 = 1; N_a1a2 = 1; N_a1a3 = 1; N_a2a3 = 1

    msds = np.array([msds_all[1]/N_a1, msds_all[3]/N_a2, msds_all[5]/N_a3, (msds_all[0] - msds_all[1])/N_a1, (msds_all[2] - msds_all[3])/N_a2, (msds_all[4] - msds_all[5])/N_a3, msds_all[0]/N_a1, msds_all[2]/N_a2, msds_all[4]/N_a3, msds_all[6]/N_a1a2, msds_all[7]/N_a1a3, msds_all[8]/N_a2a3])
    
    #outPrint('msd_self_'+cation_name,    times, msds[0])
    #outPrint('msd_self_'+anion_name,     times, msds[1])
    outPrint('msd_self_'+polyanion_name, times, msds[2])

    #outPrint('msd_dis_'+cation_name,    times, msds[3])
    #outPrint('msd_dis_'+anion_name,     times, msds[4])
    outPrint('msd_dis_'+polyanion_name, times, msds[5])

    #outPrint('msd_all_'+cation_name,    times, msds[6])
    #outPrint('msd_all_'+anion_name,     times, msds[7])
    #outPrint('msd_all_'+polyanion_name, times, msds[8])

    #outPrint('msd_cross_'+cation_name+'_'+anion_name,     times, msds[9])
    outPrint('msd_cross_'+cation_name+'_'+polyanion_name, times, msds[10])
    outPrint('msd_cross_'+anion_name+'_'+polyanion_name,  times, msds[11])

if __name__ == "__main__":
    main(trj_file, top_file, t_min, t_max, step, nt)
