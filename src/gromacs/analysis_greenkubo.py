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

import multiprocessing as mp
import functools

import MDAnalysis as mda
import numpy as np

def outPrint(filename, dt, msd):  # print msd
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# time msd', file=anaout)
        for i in range(0, len(dt)):
            print('{:10.3f} {:10.5f}'.format(dt[i], msd[i]), file=anaout)


def autocorrFFT(x):
    """ 
    Calculates the position autocorrelation function using the fast Fourier transform.
    """
    N=len(x)
    F = np.fft.fft(x, n=2*N)  
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res= (res[:N]).real   
    n=N*np.ones(N)-np.arange(0,N) 
    return res/n 


def msd_fft(r):
    """ 
    Calculates mean square displacement of the array r using the fast Fourier transform.
    """
    N=len(r)
    D=np.square(r).sum(axis=1) 
    D=np.append(D,0) 
    S2=sum([autocorrFFT(r[:, i]) for i in range(r.shape[1])])
    Q=2*D.sum()
    S1=np.zeros(N)
    for m in range(N):
        Q=Q-D[m-1]-D[N-m]
        S1[m]=Q/(N-m)
    return S1-2*S2


def calc_cond(u, anions, cations):
    """ 
    Calculates the time-averaged "mean square displaement" (the term in brackets in the Einstein relation for conductivity) given an MDAnalysis universe (u) and a selection of atoms or molecules (sel)
    
    Args: 
    u: MDAnalysis Universe
    sel: Selection of atoms (an MDAnalysis AtomGroup). Current code assumes anion and cation
    selections are single atoms, not molecules.
    
    Returns an array of "MSD" values for each time in the trajectory. 
    
    NOTE: coordinates must be unwrapped (in dcd file when creating MDAnalysis Universe)
    
    NOTE: current code assumes ions are monovalent.
    """
    qr = []
    for ts in u.trajectory:
        qr_temp = np.zeros(3)
        for anion in anions.atoms:
            qr_temp += -anion.position
        for cation in cations.atoms:
            qr_temp += cation.position
        qr.append(qr_temp)
    msd = (msd_fft(np.array(qr)))
    return msd


def cross_corr(x, y):
    N=len(x)
    F1 = np.fft.fft(x, n=2**(N*2 - 1).bit_length())  #2*N because of zero-padding, use next highest power of 2
    F2 = np.fft.fft(y, n=2**(N*2 - 1).bit_length())
    PSD = F1 * F2.conjugate()
    res = np.fft.ifft(PSD)
    res= (res[:N]).real   
    n=N*np.ones(N)-np.arange(0,N) 
    return res/n 


def msd_fft_cross(r, k):
    """ 
    Calculates mean square displacement of the array r using the fast Fourier transform.
    """
    N=len(r)
    D=np.multiply(r,k).sum(axis=1) 
    D=np.append(D,0) 
    S2=sum([cross_corr(r[:, i], k[:,i]) for i in range(r.shape[1])])
    S3=sum([cross_corr(k[:, i], r[:,i]) for i in range(k.shape[1])])
    Q=2*D.sum()
    S1=np.zeros(N)
    for m in range(N):
        Q=Q-D[m-1]-D[N-m]
        S1[m]=Q/(N-m)
    return S1-S2-S3


def calc_NE(u, atoms, times):
    """ 
    Calculates Nernst Einstein (self) conductivity for either cations or anions.
    Atoms should be the list of all cations or all anions
    Anions and cations are single atoms. 
    """
    # Current code assumes anion and cation selections are single atoms
    # Also assumes monovalent ions
    cond_NE = np.zeros(len(times))
    for atom in (atoms.atoms):
        r = []
        for ts in u.trajectory:
            r_temp = atom.position - u.atoms.center_of_mass()
            r.append(r_temp)
        msd_temp = msd_fft(np.array(r))
        cond_NE += msd_temp
    return cond_NE


def calc_catCatCorr(u, atoms, times):
    """ 
    Calculates cation-cation or anion-anion correlations.
    Atoms should be the list of all cations or anions 
    """
    r = []
    for ts in (u.trajectory):
        r_temp = np.zeros([3, atoms.n_atoms])
        for atom_num in (range((atoms.n_atoms))):
            atom = atoms[atom_num]
            r_temp[:,atom_num] += (atom.position - u.atoms.center_of_mass()) 
        r.append(r_temp)
    msd = np.zeros(len(times))
    for atom1_num in (range(atoms.n_atoms)):
        for atom2_num in (range(atoms.n_atoms)):
            if atoms[atom1_num].id != atoms[atom2_num].id:
                msd += msd_fft_cross(np.array(r)[:,:,atom1_num],np.array(r)[:,:,atom2_num])
    return msd


def calc_catAnCorr(u, cations, anions, times):
    """ 
    Calculates cation-anion correlations. Returned result is 2*(cation-anion distinct
    conductivity), accounting for the fact that sigma_(+-) = sigma_(-+).
    """
    r_cat = []
    r_an = []
    for ts in (u.trajectory):
        r_cat_temp = np.zeros(3)
        r_an_temp = np.zeros(3)
        for cation in (cations.atoms):
            r_cat_temp+=(cation.position - u.atoms.center_of_mass())  
        for anion in anions.atoms: 
            r_an_temp+=(anion.position - u.atoms.center_of_mass())
        r_cat.append(r_cat_temp)
        r_an.append(r_an_temp)
    msd = -2*msd_fft_cross(np.array(r_cat),np.array(r_an))
    return msd


def main():
    # Unit Conversions
    A2cm = 1e-8
    ps2s = 1e-12
    fs2ps = 1e-3
    e2c = 1.60217662e-19  
    
    # Constants
    q = 1.6e-19  # C, elementary charge
    kb = 1.3806504e-23  # J/K, Boltzmann constant
    T = 300
    F = 96485  # C/mol e-, Faraday's constant
    top = 'cg_topol.tpr'
    trj = 'unwrap.xtc'
    uta = mda.Universe(top, trj)
    ag1 = uta.select_atoms("type PF")
    ag2 = uta.select_atoms("type IM")
    timestep = 1 # number of femtoseconds between data collection 
    times = []
    step = 0
    run_start = 0
    for ts in uta.trajectory[:]:
        times.append(step * timestep * fs2ps)
        step += 1
    times = np.array(times)
    msd_self_ag1 = calc_NE(uta, ag1, times)
    outPrint('msd_self_ag1', times, msd_self_ag1)
    msd_self_ag2 = calc_NE(uta, ag2, times)
    outPrint('msd_self_ag2', times, msd_self_ag2)
    msd_dis_ag11 = calc_catCatCorr(uta, ag1, times)
    outPrint('msd_dis_ag11', times, msd_dis_ag11)
    msd_dis_ag22 = calc_catCatCorr(uta, ag2, times)
    outPrint('msd_dis_ag22', times, msd_dis_ag22)
    msd_dis_ag12 = calc_catAnCorr(uta, ag1, ag2, times)
    outPrint('msd_dis_ag12', times, msd_dis_ag12)

if __name__ == "__main__":
    main()