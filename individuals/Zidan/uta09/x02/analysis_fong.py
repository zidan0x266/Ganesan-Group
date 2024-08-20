"""
Copyright (C) 2018-2022 Zidan Zhang <zhangzidan@gmail.com>

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

import multiprocessing as mp
import functools

import numpy as np
import MDAnalysis as mda
from scipy import stats

def outPrint(filename, dt, msd):  # print msd
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# time msd', file=anaout)
        for i in range(0, len(dt)):
            print('{:10.3f} {:10.5f}'.format(dt[i] / 1000, msd[i]), file=anaout)


def create_mda(top, trj):
    """
    Creates MDAnalysis universe with trajectory data.
    :param top: string, path to Gromacs topology file with atom coordinates and topology
    :param trj: string or list[string], path(s) to Gromacs xtc files with trajectory data
    :return uta: MDAnalysis universe
    """
    uta = mda.Universe(top, trj)
    return uta


def define_atom_types(uta, anion_type = "type TF", cation_type = "type IM"):
    """ 
    Sorts atoms in the MDAnalysis universe based on type (cation and anion).
    Selections must be a single atom (rather than a molecule center of mass).
    :param uta: MDAnalysis universe
    :param anion_type: string, type number corresponding to anions in the Gromacs topology files
    :param cation_type: string, type number corresponding to cations in the Gromacs topology files
    :return cations, anions: MDAnalysis AtomGroups corresponding to cations and anions
    """
    anions = uta.select_atoms(anion_type)
    cations = uta.select_atoms(cation_type)
    return cations, anions


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


def create_position_arrays(uta, fun, anions, cations, times):
    """
    Creates an array containing the positions of all cations and anions over time.
    :param uta: MDAnalysis universe
    :param anions: MDAnalysis AtomGroup containing all anions (assumes anions are single atoms)
    :param cations: MDAnalysis AtomGroup containing all cations (assumes cations are single atoms)
    :param times: array[float], times at which position data was collected in the simulation
    :return anion_positions, cation_positions: array[float,float,float], array with all
    cation/anion positions. Indices correspond to time, ion index, and spatial dimension
    (x,y,z), respectively
    """
    time = 0
    gC = uta.select_atoms('all')
    #gP = uta.select_atoms('type IM or type BB or type BZ or type BS or type ST')
    gP = uta.select_atoms('type IM or type BB or type BZ')
    anion_positions = np.zeros((len(times), len(anions), 3))
    cation_positions = np.zeros((len(times), len(cations), 3))
    for ts in uta.trajectory:
        if fun == "com":
            ref = gC.center_of_mass()
        elif fun == "pom":
            ref = gP.center_of_mass()
        anion_positions[time, :, :] = anions.positions - ref
        cation_positions[time, :, :] = cations.positions - ref
        time += 1
    return anion_positions, cation_positions


def create_position_arrays_dimensional(uta, fun, anions, cations, times, axis):
    """
    Creates an array containing the positions of all cations and anions over time.
    :param uta: MDAnalysis universe
    :param anions: MDAnalysis AtomGroup containing all anions (assumes anions are single atoms)
    :param cations: MDAnalysis AtomGroup containing all cations (assumes cations are single atoms)
    :param times: array[float], times at which position data was collected in the simulation
    :return anion_positions, cation_positions: array[float,float,float], array with all
    cation/anion positions. Indices correspond to time, ion index, and spatial dimension
    (x,y,z), respectively
    """
    time = 0
    gC = uta.select_atoms('all')
    gP = uta.select_atoms('type IM or type BB or type BZ or type BS or type ST')
    #gP = uta.select_atoms('type IM or type BB or type BZ')
    anion_positions = np.zeros((len(times), len(anions), 2))
    cation_positions = np.zeros((len(times), len(cations), 2))
    for ts in uta.trajectory:
        if fun == "com":
            ref = gC.center_of_mass()
        elif fun == "pom":
            ref = gP.center_of_mass()
        tmp_anion = anions.positions - ref
        tmp_cation = cations.positions - ref
        para_anion = np.delete(tmp_anion, axis, 1)  # 0: row; 1: column
        para_cation = np.delete(tmp_cation, axis, 1)
        anion_positions[time, :, :] = para_anion
        cation_positions[time, :, :] = para_cation
        time += 1
    return anion_positions, cation_positions


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


def calc_Lij(cation_positions, anion_positions, times):
    """
    Calculates the "MSD" for the off-diagonal transport coefficient L^{ij}, i \neq j.
    :param cation_positions, anion_positions: array[float,float,float], position of each 
    atom (anion or cation, respectively) over time. Indices correspond to time, ion index,
    and spatial dimension (x,y,z), respectively.
    :param times: array[float], times at which position data was collected in the simulation
    :return msd: array[float], "MSD" corresponding to the L^{ij} transport coefficient at 
    each time.
    """
    r_cat = np.sum(cation_positions, axis = 1)
    r_an = np.sum(anion_positions, axis = 1)
    msd = msd_fft_cross(np.array(r_cat), np.array(r_an))
    return np.array(msd)

    
def compute_all_Lij(cation_positions, anion_positions, times):
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
    msd_self_cation = calc_Lii_self(cation_positions, times)
    msd_self_anion =  calc_Lii_self(anion_positions, times)
    msd_cation = calc_Lii(cation_positions, times)
    msd_anion = calc_Lii(anion_positions, times)
    msd_distinct_catan = calc_Lij(cation_positions, anion_positions, times)
    msds_all = [msd_cation, msd_self_cation, msd_anion, msd_self_anion, msd_distinct_catan]
    return msds_all


def main():
    top = 'cg_topol.tpr'
    trj = 'unwrap.xtc'
    fun = 'TEMPLATE'
    axis = 1
    uta = create_mda(top, trj)
    cations, anions = define_atom_types(uta)
    t_series = np.linspace(0, len(uta.trajectory) - 1, len(uta.trajectory))
    #anion_positions, cation_positions = create_position_arrays(uta, fun, anions, cations, t_series)
    anion_positions, cation_positions = create_position_arrays_dimensional(uta, fun, anions, cations, t_series, axis)
    msds_all = compute_all_Lij(cation_positions, anion_positions, t_series)
    outPrint('Lsi_TEMPLATE', t_series, msds_all[3])
    outPrint('Lsj_TEMPLATE', t_series, msds_all[1])
    outPrint('Lii_TEMPLATE', t_series, msds_all[2])
    outPrint('Ljj_TEMPLATE', t_series, msds_all[0])
    outPrint('Lij_TEMPLATE', t_series, msds_all[4])


if __name__ == "__main__":
    main()