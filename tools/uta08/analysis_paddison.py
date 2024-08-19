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
import numpy as np
import h5py
import math
import time
from collections import Counter
import MDAnalysis as mda


def vghopPrint(filename, shift, F1, F2, F3, F4):  # Venkat Ganesan method
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# T1, T2, T3, T4, delta_t = {}'.format(shift), file=anaout)
        for i in range(len(F1)):
            print('{:8.5f} {:8.5f} {:8.5f} {:8.5f}'.format(F1[i], F2[i], F3[i], F4[i]), file=anaout)


def sphopPrint(filename, shift, F1, F2):  # Stephen Paddison method
    with open('{}.xvg'.format(filename), 'w') as anaout:
        print('# T1, T2, delta_t = {}'.format(shift), file=anaout)
        for i in range(len(F1)):
            print('{:8.5f} {:8.5f}'.format(F1[i], F2[i]), file=anaout)


def assoFreq(frequence):
    """
    This function gets the normolized probability for association events.
    """
    fNum = frequence[0]
    for i in range(1, len(frequence)):
        fNum += frequence[i]
    numAsso = {}
    for k, v in fNum.items():
        numAsso[k] = v / float(len(frequence))
    totalAsso = sum(numAsso.values())
    avgAsso = {}
    for k, v in numAsso.items():
        avgAsso[k] = v / float(totalAsso)
    pairNum = []
    pairPro = []
    for k in sorted(avgAsso.keys()):
        pairNum.append(k)
        pairPro.append(avgAsso[k])
    return pairNum, pairPro


def anaPrint(filename, arrayNum, arrayPro):  # print association properties
    anaout = open(filename + '.xvg', 'w')
    print('# ' + filename + ' Probability', file=anaout)
    for i in range(0, len(arrayNum)):
        print('{} {:8.5f}'.format(arrayNum[i], arrayPro[i]), file=anaout)
    anaout.close()


def ds_append(ds, val):
    ds.resize(ds.shape[0]+1, axis=0)
    ds[-1] = val


def rawLoad(h5, dataset_name, mode='r'):
    if isinstance(h5, str):
        h5 = h5py.File(open(h5, 'rb'), mode)
        return h5[dataset_name], h5
    return h5[dataset_name]


def construction(h5, dataset_name, length, dtype):
    dtype = h5py.vlen_dtype(np.dtype('int32'))
    if dataset_name not in h5:
        #return h5.create_dataset(dataset_name, (0, length), dtype=dtype)
        return h5.create_dataset(dataset_name, (0, length), maxshape=(None, length), dtype=dtype, chunks=(100, length))
    return h5[dataset_name]


def get_solvation(ref_atom, ASSO, DP):
    atom = ref_atom
    sol_shell, chain_id = [], []
    for pair in range(len(ASSO)):
        atom1, atom2 = ASSO[pair]
        if atom1 == atom:
            sol_shell.append(atom2)
            chain_id.append(int(math.floor(int(atom2) / DP) + 1))
    return list(set(sol_shell)), list(set(chain_id))


def solvation(pst, Nref_atom, DP, frame_id):
    """
    get the first solvation shell for each lithium
    """
    H5 = h5py.File(pst, 'r')
    ASSO = rawLoad(H5, 'htt')
    Sol_shell, Chain_info = [], []
    for ref_atom in range(Nref_atom):
        sol_shell, chain_id = get_solvation(ref_atom, ASSO[frame_id], DP)
        Sol_shell.append(sol_shell)
        Chain_info.append(chain_id)
    return Sol_shell, Chain_info


def generate_solvation(filename, pst, Nref_atom, DP, nframes, nt):
    dtype = h5py.vlen_dtype(np.dtype('int32'))
    h5 = h5py.File(filename + '.h5', 'w')  # save solvation shell information
    Sol_shell = construction(h5, 'shell', Nref_atom, dtype)
    Chain_info = construction(h5, 'chain', Nref_atom, dtype)
    frame_ids = [frame for frame in range(nframes)]
    solvation_analyse = functools.partial(solvation, pst, Nref_atom, DP)
    pool = mp.Pool(processes=nt)
    data = pool.imap(solvation_analyse, frame_ids)
    for Shell, Chain in data:
        ds_append(Sol_shell, Shell)
        ds_append(Chain_info, Chain)
    h5.close()


def vghop(filename, h5flag1, h5flag2, shift, ref_frame):
    Vehicle, h5 = rawLoad(filename + '.h5', h5flag1)
    Chain, h5 = rawLoad(filename + '.h5', h5flag2)
    t1, t2, t3, t4 = 0, 0, 0, 0
    for i in range(len(Vehicle[ref_frame + shift])):
        if len(np.setdiff1d(Vehicle[ref_frame + shift][i], Vehicle[ref_frame][i])) == 0:
            t4 += 1
        elif len(np.setdiff1d(Vehicle[ref_frame + shift][i], Vehicle[ref_frame][i])) != 0 and len(np.setdiff1d(Chain[ref_frame + shift][i], Chain[ref_frame][i])) == 0:
            t1 += 1
        elif len(np.setdiff1d(Vehicle[ref_frame + shift][i], Vehicle[ref_frame][i])) != 0 and len(np.setdiff1d(Chain[ref_frame + shift][i], Chain[ref_frame][i])) != 0:
            t2 += 1
        elif len(Chain[ref_frame + shift][i]) == 0:
            t3 += 1
    total = t1 + t2 + t3 + t4
    return t1 / total, t2 / total, t3 / total, t4 / total


def sphop(filename, h5flag1, h5flag2, top, trj, atom_type, r_star, shift, ref_frame):
    uta = mda.Universe(top, trj)
    atg = uta.select_atoms("type " + atom_type)  # atom group
    ts = uta.trajectory[ref_frame]
    ref_position = atg.positions
    ts = uta.trajectory[ref_frame + shift]
    cur_position = atg.positions
    displ = cur_position - ref_position
    displacement = np.sqrt(np.sum(np.square(displ), axis = 1))
    mobile, immobile = [], []
    for index, length in enumerate(displacement):
        if length < r_star:
            immobile.append(index)
        else:
            mobile.append(index)
    Vehicle, h5 = rawLoad(filename + '.h5', h5flag1)
    Chain, h5 = rawLoad(filename + '.h5', h5flag2)
    t1, t2 = 0, 0
    asso_m, chain_m = [], []
    asso_t, chain_t = [], []  # total
    asso_im, chain_im = [], []
    for i in range(len(Vehicle[ref_frame + shift])):
        if i in mobile:
            if len(np.setdiff1d(Vehicle[ref_frame + shift][i], Vehicle[ref_frame][i])) != 0 and len(np.setdiff1d(Chain[ref_frame + shift][i], Chain[ref_frame][i])) == 0:
                t1 += 1
            elif len(np.setdiff1d(Vehicle[ref_frame + shift][i], Vehicle[ref_frame][i])) != 0 and len(np.setdiff1d(Chain[ref_frame + shift][i], Chain[ref_frame][i])) != 0:
                t2 += 1
            asso_m.append(len(Vehicle[ref_frame][i]))
            chain_m.append(len(Chain[ref_frame][i]))
            asso_t.append(len(Vehicle[ref_frame][i]))
            chain_t.append(len(Chain[ref_frame][i]))
        else:
            asso_im.append(len(Vehicle[ref_frame][i]))
            chain_im.append(len(Chain[ref_frame][i]))
            asso_t.append(len(Vehicle[ref_frame][i]))
            chain_t.append(len(Chain[ref_frame][i]))
    total = t1 + t2
    return t1 / total, t2 / total, Counter(asso_m), Counter(chain_m), Counter(asso_im), Counter(chain_im), Counter(asso_t), Counter(chain_t)


def main():
    # generic parameters
    func = 1  # if 0: generate the .h5 file, else using exsiting .h5 file
    hop_style = 'sp_parallel'
    nframes = 200001  # number of frames
    npoints = 1001
    nt = 56  # number of processors
    DP = 32
    Natoms = 256
    shift = 1
    pst = 'pils.h5'
    top = '../cg_topol.tpr'
    trj = '../fong/unwrap.xtc'
    atom_type = 'PF'                           # fetch central particle
    h5file = 'solvation'
    h5flag1 = 'shell'
    h5flag2 = 'chain'
    r_star = 4.0  # \AA obtained from the van Hove function
    if func == 0:
        generate_solvation(h5file, pst, Natoms, DP, nframes, nt)
    if hop_style == 'vg':
        F1, F2, F3, F4 = [], [], [], []
        for frame in np.linspace(0, nframes, npoints, dtype = int):
            if frame + shift <= nframes - 1:
                hop_types = vghop(h5file, h5flag1, h5flag2, shift, frame)
                F1.append(hop_types[0])
                F2.append(hop_types[1])
                F3.append(hop_types[2])
                F4.append(hop_types[3])
        vghopPrint('vg_hopping_{}'.format(shift), shift, F1, F2, F3, F4)
    elif hop_style == 'vg_parallel':
        frame_ids = [frame for frame in np.linspace(0, nframes, npoints, dtype = int) if (frame + shift) <= (nframes - 1)]
        vghop_analysis = functools.partial(vghop, h5file, h5flag1, h5flag2, shift)
        pool = mp.Pool(processes=nt)
        data = pool.imap(vghop_analysis, frame_ids)
        F1, F2, F3, F4 = [], [], [], []
        for f1, f2, f3, f4 in data:
            F1.append(f1)
            F2.append(f2)
            F3.append(f3)
            F4.append(f4)
        vghopPrint('vg_hopping_{}'.format(shift), shift, F1, F2, F3, F4)
    elif hop_style == 'sp':
        F1, F2 = [], []
        freq_asso_m, freq_chain_m = [], []
        freq_asso_t, freq_chain_t = [], []
        freq_asso_im, freq_chain_im = [], []
        for frame in np.linspace(0, nframes, npoints, dtype = int):
            if frame + shift <= nframes:
                hop_types = sphop(h5file, h5flag1, h5flag2, top, trj, atom_type, r_star, shift, frame)
                F1.append(hop_types[0])
                F2.append(hop_types[1])
                freq_asso_m.append(hop_types[2])
                freq_chain_m.append(hop_types[3])
                freq_asso_im.append(hop_types[4])
                freq_chain_im.append(hop_types[5])
                freq_asso_t.append(hop_types[6])
                freq_chain_t.append(hop_types[7])
        sphopPrint('sp_hopping_{}'.format(shift), shift, F1, F2)
        asso_Num_m, asso_Pro_m = assoFreq(freq_asso_m)
        anaPrint('assoPair_m', asso_Num_m, asso_Pro_m)
        chain_Num_m, chain_Pro_m = assoFreq(freq_chain_m)
        anaPrint('assoChain_m', chain_Num_m, chain_Pro_m)
        asso_Num_im, asso_Pro_im = assoFreq(freq_asso_im)
        anaPrint('assoPair_im', asso_Num_im, asso_Pro_im)
        chain_Num_im, chain_Pro_im = assoFreq(freq_chain_im)
        anaPrint('assoChain_im', chain_Num_im, chain_Pro_im)
        asso_Num_t, asso_Pro_t = assoFreq(freq_asso_t)
        anaPrint('assoPair_t', asso_Num_t, asso_Pro_t)
        chain_Num_t, chain_Pro_t = assoFreq(freq_chain_t)
        anaPrint('assoChain_t', chain_Num_t, chain_Pro_t)
    elif hop_style == 'sp_parallel':
        frame_ids = [frame for frame in np.linspace(0, nframes, npoints, dtype = int) if (frame + shift) <= (nframes - 1)]
        sphop_analysis = functools.partial(sphop, h5file, h5flag1, h5flag2, top, trj, atom_type, r_star, shift)
        pool = mp.Pool(processes=nt)
        data = pool.imap(sphop_analysis, frame_ids)
        F1, F2 = [], []
        freq_asso_m, freq_chain_m = [], []
        freq_asso_t, freq_chain_t = [], []
        freq_asso_im, freq_chain_im = [], []
        for f1, f2, asso_m, chain_m, asso_im, chain_im, asso_t, chain_t in data:
            F1.append(f1)
            F2.append(f2)
            freq_asso_m.append(asso_m)
            freq_chain_m.append(chain_m)
            freq_asso_im.append(asso_im)
            freq_chain_im.append(chain_im)
            freq_asso_t.append(asso_t)
            freq_chain_t.append(chain_t)
        sphopPrint('sp_hopping_{}'.format(shift), shift, F1, F2)
        asso_Num_m, asso_Pro_m = assoFreq(freq_asso_m)
        anaPrint('assoPair_m', asso_Num_m, asso_Pro_m)
        chain_Num_m, chain_Pro_m = assoFreq(freq_chain_m)
        anaPrint('assoChain_m', chain_Num_m, chain_Pro_m)
        asso_Num_im, asso_Pro_im = assoFreq(freq_asso_im)
        anaPrint('assoPair_im', asso_Num_im, asso_Pro_im)
        chain_Num_im, chain_Pro_im = assoFreq(freq_chain_im)
        anaPrint('assoChain_im', chain_Num_im, chain_Pro_im)
        asso_Num_t, asso_Pro_t = assoFreq(freq_asso_t)
        anaPrint('assoPair_t', asso_Num_t, asso_Pro_t)
        chain_Num_t, chain_Pro_t = assoFreq(freq_chain_t)
        anaPrint('assoChain_t', chain_Num_t, chain_Pro_t)
        

if __name__ == "__main__":
    main()
