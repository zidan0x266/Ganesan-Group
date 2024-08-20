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

import math
import MDAnalysis as mda
import MDAnalysis.analysis.distances as dist
import MDAnalysis.lib.distances as distances
import numpy as np
from collections import Counter


def ioninDPrint(filename, ioninD, withsalt):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        #print('# data explanation', file=anaout)
        if withsalt:
            for i in range(len(ioninD)):
                print('{:5d} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d} {:5d}'.format(ioninD[i][0], ioninD[i][1], ioninD[i][2], ioninD[i][3], ioninD[i][4], ioninD[i][5], ioninD[i][6], ioninD[i][7], ioninD[i][8]), file=anaout)
        else:
            for i in range(len(ioninD)):
                print('{:5d} {:5d} {:5d} {:5d} {:5d} {:5d}'.format(ioninD[i][0], ioninD[i][1], ioninD[i][2], ioninD[i][3], ioninD[i][4], ioninD[i][5]), file=anaout)


def assoFreq(frequency):
    """
    This function gets the normalized probability for association events.
    """
    reservoir = frequency[0]
    for i in range(1, len(frequency)):
        reservoir += frequency[i]
    counts = {}
    for k, v in reservoir.items():
        counts[k] = v / float(len(frequency))  # average over samples
    total = sum(counts.values())
    prob = {}
    for k, v in counts.items():
        prob[k] = v / float(total)  # normalization
    Num = []
    Pro = []
    for k in sorted(prob.keys()):
        Num.append(k)
        Pro.append(prob[k])
    return Num, Pro


def assoPrint(filename, arrayNum, arrayPro):  # print association properties
    anaout = open(filename + '.xvg', 'w')
    print('# ' + filename + ' Probability', file=anaout)
    for i in range(0, len(arrayNum)):
        print('{} {:10.5f}'.format(arrayNum[i], arrayPro[i]), file=anaout)
    anaout.close()


def get_domains(concentration):
    if concentration == 'c000':
        domains = [29.387, 9.635, 97.136, 10.542]  # c000
    elif concentration == 'c020':
        domains = [28.330, 11.155, 100.320, 10.950]  # c020
    elif concentration == 'c040':
        domains = [26.760, 10.590, 102.400, 11.142]  # c040
    elif concentration == 'c060':
        domains = [26.140, 10.776, 105.340, 10.446]  # c060
    elif concentration == 'c080':
        domains = [25.320, 10.667, 107.740, 11.429]  # c080
    # define regions
    if1 = [domains[0] - domains[1] / 2, domains[0] + domains[1] / 2]
    bk1 = [domains[0] + domains[1] / 2, domains[2] - domains[3] / 2]
    if2 = [domains[2] - domains[3] / 2, domains[2] + domains[3] / 2]
    to1 = [domains[0] - domains[1] / 2, domains[2] + domains[3] / 2]
    return if1, bk1, if2, to1


def get_sub_domains(concentration):
    if concentration == 'c000':
        domains = [29.387, 9.635, 97.136, 10.542]  # c000
    elif concentration == 'c020':
        domains = [28.330, 11.155, 100.320, 10.950]  # c020
    elif concentration == 'c040':
        domains = [26.760, 10.590, 102.400, 11.142]  # c040
    elif concentration == 'c060':
        domains = [26.140, 10.776, 105.340, 10.446]  # c060
    elif concentration == 'c080':
        domains = [25.320, 10.667, 107.740, 11.429]  # c080
    # define sub-regions
    bk1 = [domains[0] + domains[1] / 2, domains[2] - domains[3] / 2]
    length_bk = bk1[1] - bk1[0]
    increment = float("{0:.3f}".format(length_bk / 3))  # split the bulk region into three sub-regions
    sbkl = [bk1[0], bk1[0] + increment]
    sbkm = [bk1[0] + increment, bk1[0] + 2 * increment]
    sbkr = [bk1[0] + 2 * increment, bk1[1]]
    return sbkl, sbkm, sbkr


def get_ions_inD(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
    if withsalt:
        group3 = uta.select_atoms("type " + ion3)
    ts = uta.trajectory[frame_id]
    if1, bk1, if2, to1 = get_domains(concentration)
    ag1inif1 = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= if1[0] and group1.positions[:, dimension][i] < if1[1])]
    ag1inif2 = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= if2[0] and group1.positions[:, dimension][i] <= if2[1])]
    ag1inbk = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= bk1[0] and group1.positions[:, dimension][i] < bk1[1])]
    ag1into = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= to1[0] and group1.positions[:, dimension][i] <= to1[1])]
    ag2inif1 = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= if1[0] and group2.positions[:, dimension][i] < if1[1])]
    ag2inif2 = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= if2[0] and group2.positions[:, dimension][i] <= if2[1])]
    ag2inbk = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= bk1[0] and group2.positions[:, dimension][i] < bk1[1])]
    ag2into = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= to1[0] and group2.positions[:, dimension][i] <= to1[1])]
    if withsalt:
        ag3inif1 = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= if1[0] and group3.positions[:, dimension][i] < if1[1])]
        ag3inif2 = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= if2[0] and group3.positions[:, dimension][i] <= if2[1])]
        ag3inbk = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= bk1[0] and group3.positions[:, dimension][i] < bk1[1])]
        ag3into = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= to1[0] and group3.positions[:, dimension][i] <= to1[1])]
    ag1inif = ag1inif1 + ag1inif2
    ag2inif = ag2inif1 + ag2inif2
    if withsalt:
        ag3inif = ag3inif1 + ag3inif2
        data = [len(ag1into), len(ag1inbk), len(ag1inif), len(ag2into), len(ag2inbk), len(ag2inif), len(ag3into), len(ag3inbk), len(ag3inif)]
    else:
        data = [len(ag1into), len(ag1inbk), len(ag1inif), len(ag2into), len(ag2inbk), len(ag2inif)]
    return ts.frame, data


def get_ions_insubD(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
    if withsalt:
        group3 = uta.select_atoms("type " + ion3)
    ts = uta.trajectory[frame_id]
    sbkl, sbkm, sbkr = get_sub_domains(concentration)
    ag1inr1 = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= sbkl[0] and group1.positions[:, dimension][i] < sbkl[1])]
    ag1inr2 = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= sbkm[0] and group1.positions[:, dimension][i] < sbkm[1])]
    ag1inr3 = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= sbkr[0] and group1.positions[:, dimension][i] <= sbkr[1])]
    ag2inr1 = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= sbkl[0] and group2.positions[:, dimension][i] < sbkl[1])]
    ag2inr2 = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= sbkm[0] and group2.positions[:, dimension][i] < sbkm[1])]
    ag2inr3 = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= sbkr[0] and group2.positions[:, dimension][i] <= sbkr[1])]
    if withsalt:
        ag3inr1 = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= sbkl[0] and group3.positions[:, dimension][i] < sbkl[1])]
        ag3inr2 = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= sbkm[0] and group3.positions[:, dimension][i] < sbkm[1])]
        ag3inr3 = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= sbkr[0] and group3.positions[:, dimension][i] <= sbkr[1])]
    if withsalt:
        data = [len(ag1inr1), len(ag1inr2), len(ag1inr3), len(ag2inr1), len(ag2inr2), len(ag2inr3), len(ag3inr1), len(ag3inr2), len(ag3inr3)]
    else:
        data = [len(ag1inr1), len(ag1inr2), len(ag1inr3), len(ag2inr1), len(ag2inr2), len(ag2inr3)]
    return ts.frame, data


def get_ionlist_inD(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
    if withsalt:
        group3 = uta.select_atoms("type " + ion3)
    ts = uta.trajectory[frame_id]
    if1, bk1, if2, to1 = get_domains(concentration)
    ag1inif1 = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= if1[0] and group1.positions[:, dimension][i] < if1[1])]
    ag1inif2 = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= if2[0] and group1.positions[:, dimension][i] <= if2[1])]
    ag1inbk = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= bk1[0] and group1.positions[:, dimension][i] < bk1[1])]
    ag1into = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= to1[0] and group1.positions[:, dimension][i] <= to1[1])]
    ag2inif1 = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= if1[0] and group2.positions[:, dimension][i] < if1[1])]
    ag2inif2 = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= if2[0] and group2.positions[:, dimension][i] <= if2[1])]
    ag2inbk = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= bk1[0] and group2.positions[:, dimension][i] < bk1[1])]
    ag2into = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= to1[0] and group2.positions[:, dimension][i] <= to1[1])]
    if withsalt:
        ag3inif1 = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= if1[0] and group3.positions[:, dimension][i] < if1[1])]
        ag3inif2 = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= if2[0] and group3.positions[:, dimension][i] <= if2[1])]
        ag3inbk = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= bk1[0] and group3.positions[:, dimension][i] < bk1[1])]
        ag3into = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= to1[0] and group3.positions[:, dimension][i] <= to1[1])]
    ag1inif = ag1inif1 + ag1inif2
    ag2inif = ag2inif1 + ag2inif2
    if withsalt:
        ag3inif = ag3inif1 + ag3inif2
        return ag1into, ag1inbk, ag1inif, ag2into, ag2inbk, ag2inif, ag3into, ag3inbk, ag3inif
    else:
        return ag1into, ag1inbk, ag1inif, ag2into, ag2inbk, ag2inif


def get_ionlist_insubD(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
    if withsalt:
        group3 = uta.select_atoms("type " + ion3)
    ts = uta.trajectory[frame_id]
    sbkl, sbkm, sbkr = get_sub_domains(concentration)
    ag1inr1 = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= sbkl[0] and group1.positions[:, dimension][i] < sbkl[1])]
    ag1inr2 = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= sbkm[0] and group1.positions[:, dimension][i] < sbkm[1])]
    ag1inr3 = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= sbkr[0] and group1.positions[:, dimension][i] <= sbkr[1])]
    ag2inr1 = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= sbkl[0] and group2.positions[:, dimension][i] < sbkl[1])]
    ag2inr2 = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= sbkm[0] and group2.positions[:, dimension][i] < sbkm[1])]
    ag2inr3 = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= sbkr[0] and group2.positions[:, dimension][i] <= sbkr[1])]
    if withsalt:
        ag3inr1 = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= sbkl[0] and group3.positions[:, dimension][i] < sbkl[1])]
        ag3inr2 = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= sbkm[0] and group3.positions[:, dimension][i] < sbkm[1])]
        ag3inr3 = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= sbkr[0] and group3.positions[:, dimension][i] <= sbkr[1])]
    if withsalt:
        return ag1inr1, ag1inr2, ag1inr3, ag2inr1, ag2inr2, ag2inr3, ag3inr1, ag3inr2, ag3inr3
    else:
        return ag1inr1, ag1inr2, ag1inr3, ag2inr1, ag2inr2, ag2inr3


def get_asso_domains_12(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, cutoff1, cutoff2, dp, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
    if withsalt:
        group3 = uta.select_atoms("type " + ion3)
        _, ag1inbk, ag1inif, _, ag2inbk, ag2inif, _, ag3inbk, ag3inif = get_ionlist_inD(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, frame_id)
    else:
        _, ag1inbk, ag1inif, _, ag2inbk, ag2inif = get_ionlist_inD(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, frame_id)
    ts = uta.trajectory[frame_id]
    cell = ts.dimensions
    asso_pair_bk, asso_chain_bk = [], []
    asso_pair_if, asso_chain_if = [], []
    for atom in ag1inbk:
        chain_ids = []
        distpair = distances.capped_distance(group1.positions[atom], group2.positions, cutoff1, box = cell)[0]
        asso_pair_bk.append(len(distpair))
        for i in range(len(distpair)):
            atom_pair = distpair[i][1]
            chain_ids.append(int(math.floor(atom_pair / dp) + 1))
        asso_chain_bk.append(len(np.unique(chain_ids)))
    for atom in ag1inif:
        chain_ids = []
        distpair = distances.capped_distance(group1.positions[atom], group2.positions, cutoff1, box = cell)[0]
        asso_pair_if.append(len(distpair))
        for i in range(len(distpair)):
            atom_pair = distpair[i][1]
            chain_ids.append(int(math.floor(atom_pair / dp) + 1))
        asso_chain_if.append(len(np.unique(chain_ids)))
    asso_pair_to = asso_pair_bk + asso_pair_if
    asso_chain_to = asso_chain_bk + asso_chain_if
    return Counter(asso_pair_to), Counter(asso_chain_to), Counter(asso_pair_bk), Counter(asso_chain_bk), Counter(asso_pair_if), Counter(asso_chain_if)


def get_asso_sub_domains_12(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, cutoff1, cutoff2, dp, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
    if withsalt:
        group3 = uta.select_atoms("type " + ion3)
        ag1inr1, ag1inr2, ag1inr3, _, _, _, ag3inr1, ag3inr2, ag3inr3 = get_ionlist_insubD(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, frame_id)
    else:
        ag1inr1, ag1inr2, ag1inr3, _, _, _ = get_ionlist_insubD(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, frame_id)
    ts = uta.trajectory[frame_id]
    cell = ts.dimensions
    asso_pair_r1, asso_chain_r1 = [], []
    asso_pair_r2, asso_chain_r2 = [], []
    asso_pair_r3, asso_chain_r3 = [], []
    for atom in ag1inr1:
        chain_ids = []
        distpair = distances.capped_distance(group1.positions[atom], group2.positions, cutoff1, box = cell)[0]
        asso_pair_r1.append(len(distpair))
        for i in range(len(distpair)):
            atom_pair = distpair[i][1]
            chain_ids.append(int(math.floor(atom_pair / dp) + 1))
        asso_chain_r1.append(len(np.unique(chain_ids)))
    for atom in ag1inr2:
        chain_ids = []
        distpair = distances.capped_distance(group1.positions[atom], group2.positions, cutoff1, box = cell)[0]
        asso_pair_r2.append(len(distpair))
        for i in range(len(distpair)):
            atom_pair = distpair[i][1]
            chain_ids.append(int(math.floor(atom_pair / dp) + 1))
        asso_chain_r2.append(len(np.unique(chain_ids)))
    for atom in ag1inr3:
        chain_ids = []
        distpair = distances.capped_distance(group1.positions[atom], group2.positions, cutoff1, box = cell)[0]
        asso_pair_r3.append(len(distpair))
        for i in range(len(distpair)):
            atom_pair = distpair[i][1]
            chain_ids.append(int(math.floor(atom_pair / dp) + 1))
        asso_chain_r3.append(len(np.unique(chain_ids)))
    return Counter(asso_pair_r1), Counter(asso_chain_r1), Counter(asso_pair_r2), Counter(asso_chain_r2), Counter(asso_pair_r3), Counter(asso_chain_r3)


def get_asso_domains_13(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, cutoff1, cutoff2, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
    if withsalt:
        group3 = uta.select_atoms("type " + ion3)
        _, ag1inbk, ag1inif, _, ag2inbk, ag2inif, _, ag3inbk, ag3inif = get_ionlist_inD(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, frame_id)
    else:
        _, ag1inbk, ag1inif, _, ag2inbk, ag2inif = get_ionlist_inD(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, frame_id)
    ts = uta.trajectory[frame_id]
    cell = ts.dimensions
    asso_pair_bk_13, asso_pair_if_13 = [], []
    asso_pair_bk_31, asso_pair_if_31 = [], []
    for atom in ag1inbk:
        distpair = distances.capped_distance(group1.positions[atom], group3.positions, cutoff2, box = cell)[0]
        asso_pair_bk_13.append(len(distpair))
    for atom in ag1inif:
        distpair = distances.capped_distance(group1.positions[atom], group3.positions, cutoff2, box = cell)[0]
        asso_pair_if_13.append(len(distpair))
    asso_pair_to_13 = asso_pair_bk_13 + asso_pair_if_13
    for atom in ag3inbk:
        distpair = distances.capped_distance(group3.positions[atom], group1.positions, cutoff2, box = cell)[0]
        asso_pair_bk_31.append(len(distpair))
    for atom in ag3inif:
        distpair = distances.capped_distance(group3.positions[atom], group1.positions, cutoff2, box = cell)[0]
        asso_pair_if_31.append(len(distpair))
    asso_pair_to_31 = asso_pair_bk_31 + asso_pair_if_31
    return Counter(asso_pair_to_13), Counter(asso_pair_bk_13), Counter(asso_pair_if_13), Counter(asso_pair_to_31), Counter(asso_pair_bk_31), Counter(asso_pair_if_31)


def get_asso_sub_domains_13(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, cutoff1, cutoff2, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
    if withsalt:
        group3 = uta.select_atoms("type " + ion3)
        ag1inr1, ag1inr2, ag1inr3, _, _, _, ag3inr1, ag3inr2, ag3inr3 = get_ionlist_inD(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, frame_id)
    else:
        ag1inr1, ag1inr2, ag1inr3, _, _, _ = get_ionlist_inD(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, frame_id)
    ts = uta.trajectory[frame_id]
    cell = ts.dimensions
    asso_pair_r1_13, asso_pair_r2_13, asso_pair_r3_13 = [], [], []
    asso_pair_r1_31, asso_pair_r2_31, asso_pair_r3_31 = [], [], []
    for atom in ag1inr1:
        distpair = distances.capped_distance(group1.positions[atom], group3.positions, cutoff2, box = cell)[0]
        asso_pair_r1_13.append(len(distpair))
    for atom in ag1inr2:
        distpair = distances.capped_distance(group1.positions[atom], group3.positions, cutoff2, box = cell)[0]
        asso_pair_r2_13.append(len(distpair))
    for atom in ag1inr3:
        distpair = distances.capped_distance(group1.positions[atom], group3.positions, cutoff2, box = cell)[0]
        asso_pair_r3_13.append(len(distpair))
    for atom in ag3inr1:
        distpair = distances.capped_distance(group3.positions[atom], group1.positions, cutoff2, box = cell)[0]
        asso_pair_r1_31.append(len(distpair))
    for atom in ag3inr2:
        distpair = distances.capped_distance(group3.positions[atom], group1.positions, cutoff2, box = cell)[0]
        asso_pair_r2_31.append(len(distpair))
    for atom in ag3inr3:
        distpair = distances.capped_distance(group3.positions[atom], group1.positions, cutoff2, box = cell)[0]
        asso_pair_r3_31.append(len(distpair))
    return Counter(asso_pair_r1_13), Counter(asso_pair_r2_13), Counter(asso_pair_r3_13), Counter(asso_pair_r1_31), Counter(asso_pair_r2_31), Counter(asso_pair_r3_31)


def main():
    missions_list = ["ions_in_domains", "ions_in_sub_domains", "ion_associations_12", "ion_associations_13", "ion_associations_sub_12", "ion_associations_sub_13"]
    missions = [True, True, True, False, True, False]
    top = "../cg_topol.tpr"
    trj = "../data_101.xtc"
    nt=20
    atom_group1 = "TF"
    atom_group2 = "IM"
    atom_group3 = "LI"
    withsalt = True
    if withsalt:
        missions[3] = True
        missions[5] = True
    dimension = 1  # 0: X; 1: Y; 2: Z
    concentration = "c060"
    cutoff1 = 7.7  # cutoff for anion - cation
    cutoff2 = 5.1  # cutoff for anion - lithium
    dp = 20  # conductive chain length
    uta = mda.Universe(top, trj)
    frame_ids = [ts.frame for ts in uta.trajectory]
    if missions[0]:
        ioninD = []
        do_analyse = functools.partial(get_ions_inD, top, trj, atom_group1, atom_group2, atom_group3, withsalt, dimension, concentration)
        pool = mp.Pool(processes=nt)
        data = pool.imap(do_analyse, frame_ids)
        for _, ion_distribution in data:
            ioninD.append(ion_distribution)
        ioninDPrint(missions_list[0], ioninD, withsalt)
        print("Analysis of ions in domains complete!")
    if missions[1]:
        ioninD = []
        do_analyse = functools.partial(get_ions_insubD, top, trj, atom_group1, atom_group2, atom_group3, withsalt, dimension, concentration)
        pool = mp.Pool(processes=nt)
        data = pool.imap(do_analyse, frame_ids)
        for _, ion_distribution in data:
            ioninD.append(ion_distribution)
        ioninDPrint(missions_list[1], ioninD, withsalt)
        print("Analysis of ions in sub-domains complete!")
    if missions[2]:
        freq_pair_to, freq_chain_to = [], []
        freq_pair_bk, freq_chain_bk = [], []
        freq_pair_if, freq_chain_if = [], []
        do_analyse = functools.partial(get_asso_domains_12, top, trj, atom_group1, atom_group2, atom_group3, withsalt, dimension, concentration, cutoff1, cutoff2, dp)
        pool = mp.Pool(processes=nt)
        data = pool.imap(do_analyse, frame_ids)
        for d1, d2, d3, d4, d5, d6 in data:
            freq_pair_to.append(d1)
            freq_chain_to.append(d2)
            freq_pair_bk.append(d3)
            freq_chain_bk.append(d4)
            freq_pair_if.append(d5)
            freq_chain_if.append(d6)
        pair_Num_to, pair_Pro_to = assoFreq(freq_pair_to)
        chain_Num_to, chain_Pro_to = assoFreq(freq_chain_to)
        assoPrint('anion_cation' + '_Pair_to', pair_Num_to, pair_Pro_to)
        assoPrint('anion_cation' + '_Chain_to', chain_Num_to, chain_Pro_to)
        pair_Num_bk, pair_Pro_bk = assoFreq(freq_pair_bk)
        chain_Num_bk, chain_Pro_bk = assoFreq(freq_chain_bk)
        assoPrint('anion_cation' + '_Pair_bk', pair_Num_bk, pair_Pro_bk)
        assoPrint('anion_cation' + '_Chain_bk', chain_Num_bk, chain_Pro_bk)
        pair_Num_if, pair_Pro_if = assoFreq(freq_pair_if)
        chain_Num_if, chain_Pro_if = assoFreq(freq_chain_if)
        assoPrint('anion_cation' + '_Pair_if', pair_Num_if, pair_Pro_if)
        assoPrint('anion_cation' + '_Chain_if', chain_Num_if, chain_Pro_if)
        print("Analysis of anion - cation associations complete!")
    if missions[3]:
        freq_pair_to_13, freq_pair_to_31 = [], []
        freq_pair_bk_13, freq_pair_bk_31 = [], []
        freq_pair_if_13, freq_pair_if_31 = [], []
        do_analyse = functools.partial(get_asso_domains_13, top, trj, atom_group1, atom_group2, atom_group3, withsalt, dimension, concentration, cutoff1, cutoff2)
        pool = mp.Pool(processes=nt)
        data = pool.imap(do_analyse, frame_ids)
        for d1, d2, d3, d4, d5, d6 in data:
            freq_pair_to_13.append(d1)
            freq_pair_bk_13.append(d2)
            freq_pair_if_13.append(d3)
            freq_pair_to_31.append(d4)
            freq_pair_bk_31.append(d5)
            freq_pair_if_31.append(d6)
        pair_Num_to_13, pair_Pro_to_13 = assoFreq(freq_pair_to_13)
        pair_Num_bk_13, pair_Pro_bk_13 = assoFreq(freq_pair_bk_13)
        pair_Num_if_13, pair_Pro_if_13 = assoFreq(freq_pair_if_13)
        assoPrint('anion_lithium' + '_Pair_to', pair_Num_to_13, pair_Pro_to_13)
        assoPrint('anion_lithium' + '_Pair_bk', pair_Num_bk_13, pair_Pro_bk_13)
        assoPrint('anion_lithium' + '_Pair_if', pair_Num_if_13, pair_Pro_if_13)
        pair_Num_to_31, pair_Pro_to_31 = assoFreq(freq_pair_to_31)
        pair_Num_bk_31, pair_Pro_bk_31 = assoFreq(freq_pair_bk_31)
        pair_Num_if_31, pair_Pro_if_31 = assoFreq(freq_pair_if_31)
        assoPrint('lithium_anion' + '_Pair_to', pair_Num_to_31, pair_Pro_to_31)
        assoPrint('lithium_anion' + '_Pair_bk', pair_Num_bk_31, pair_Pro_bk_31)
        assoPrint('lithium_anion' + '_Pair_if', pair_Num_if_31, pair_Pro_if_31)
        print("Analysis of anion - lithium associations complete!")
    if missions[4]:
        freq_pair_r1, freq_chain_r1 = [], []
        freq_pair_r2, freq_chain_r2 = [], []
        freq_pair_r3, freq_chain_r3 = [], []
        do_analyse = functools.partial(get_asso_sub_domains_12, top, trj, atom_group1, atom_group2, atom_group3, withsalt, dimension, concentration, cutoff1, cutoff2, dp)
        pool = mp.Pool(processes=nt)
        data = pool.imap(do_analyse, frame_ids)
        for d1, d2, d3, d4, d5, d6 in data:
            freq_pair_r1.append(d1)
            freq_chain_r1.append(d2)
            freq_pair_r2.append(d3)
            freq_chain_r2.append(d4)
            freq_pair_r3.append(d5)
            freq_chain_r3.append(d6)
        pair_Num_r1, pair_Pro_r1 = assoFreq(freq_pair_r1)
        chain_Num_r1, chain_Pro_r1 = assoFreq(freq_chain_r1)
        assoPrint('anion_cation' + '_Pair_r1', pair_Num_r1, pair_Pro_r1)
        assoPrint('anion_cation' + '_Chain_r1', chain_Num_r1, chain_Pro_r1)
        pair_Num_r2, pair_Pro_r2 = assoFreq(freq_pair_r2)
        chain_Num_r2, chain_Pro_r2 = assoFreq(freq_chain_r2)
        assoPrint('anion_cation' + '_Pair_r2', pair_Num_r2, pair_Pro_r2)
        assoPrint('anion_cation' + '_Chain_r2', chain_Num_r2, chain_Pro_r2)
        pair_Num_r3, pair_Pro_r3 = assoFreq(freq_pair_r3)
        chain_Num_r3, chain_Pro_r3 = assoFreq(freq_chain_r3)
        assoPrint('anion_cation' + '_Pair_r3', pair_Num_r3, pair_Pro_r3)
        assoPrint('anion_cation' + '_Chain_r3', chain_Num_r3, chain_Pro_r3)
        print("Analysis of anion - cation associations in sub-domains complete!")
    if missions[5]:
        freq_pair_r1_13, freq_pair_r1_31 = [], []
        freq_pair_r2_13, freq_pair_r2_31 = [], []
        freq_pair_r3_13, freq_pair_r3_31 = [], []
        do_analyse = functools.partial(get_asso_sub_domains_13, top, trj, atom_group1, atom_group2, atom_group3, withsalt, dimension, concentration, cutoff1, cutoff2)
        pool = mp.Pool(processes=nt)
        data = pool.imap(do_analyse, frame_ids)
        for d1, d2, d3, d4, d5, d6 in data:
            freq_pair_r1_13.append(d1)
            freq_pair_r2_13.append(d2)
            freq_pair_r3_13.append(d3)
            freq_pair_r1_31.append(d4)
            freq_pair_r2_31.append(d5)
            freq_pair_r3_31.append(d6)
        pair_Num_r1_13, pair_Pro_r1_13 = assoFreq(freq_pair_r1_13)
        pair_Num_r2_13, pair_Pro_r2_13 = assoFreq(freq_pair_r2_13)
        pair_Num_r3_13, pair_Pro_r3_13 = assoFreq(freq_pair_r3_13)
        assoPrint('anion_lithium' + '_Pair_r1', pair_Num_r1_13, pair_Pro_r1_13)
        assoPrint('anion_lithium' + '_Pair_r2', pair_Num_r2_13, pair_Pro_r2_13)
        assoPrint('anion_lithium' + '_Pair_r3', pair_Num_r3_13, pair_Pro_r3_13)
        pair_Num_r1_31, pair_Pro_r1_31 = assoFreq(freq_pair_r1_31)
        pair_Num_r2_31, pair_Pro_r2_31 = assoFreq(freq_pair_r2_31)
        pair_Num_r3_31, pair_Pro_r3_31 = assoFreq(freq_pair_r3_31)
        assoPrint('lithium_anion' + '_Pair_r1', pair_Num_r1_31, pair_Pro_r1_31)
        assoPrint('lithium_anion' + '_Pair_r2', pair_Num_r2_31, pair_Pro_r2_31)
        assoPrint('lithium_anion' + '_Pair_r3', pair_Num_r3_31, pair_Pro_r3_31)
        print("Analysis of anion - lithium associations in sub-domains complete!")


if __name__ == "__main__":
    main()