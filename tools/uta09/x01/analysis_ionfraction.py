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
from itertools import combinations


def bridgePrint1(filename, d1, d2, d3, d4):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        for i in range(len(d1)):
            print('{:8d} {:8d} {:8d} {:8d}'.format(d1[i], d2[i], d3[i], d4[i]), file=anaout)


def bridgePrint2(filename, d1):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        for i in range(len(d1)):
            print('{:8d}'.format(d1[i]), file=anaout)


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


def get_ionlist_inD(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
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
    ag1inif = ag1inif1 + ag1inif2
    ag2inif = ag2inif1 + ag2inif2
    if withsalt:
        group3 = uta.select_atoms("type " + ion3)
        ag3inif1 = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= if1[0] and group3.positions[:, dimension][i] < if1[1])]
        ag3inif2 = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= if2[0] and group3.positions[:, dimension][i] <= if2[1])]
        ag3inbk = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= bk1[0] and group3.positions[:, dimension][i] < bk1[1])]
        ag3into = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= to1[0] and group3.positions[:, dimension][i] <= to1[1])]
        ag3inif = ag3inif1 + ag3inif2
        return ag1into, ag1inbk, ag1inif, ag2into, ag2inbk, ag2inif, ag3into, ag3inbk, ag3inif
    else:
        return ag1into, ag1inbk, ag1inif, ag2into, ag2inbk, ag2inif


def get_ionlist_insubD(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
    ts = uta.trajectory[frame_id]
    sbkl, sbkm, sbkr = get_sub_domains(concentration)
    ag1inr1 = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= sbkl[0] and group1.positions[:, dimension][i] < sbkl[1])]
    ag1inr2 = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= sbkm[0] and group1.positions[:, dimension][i] < sbkm[1])]
    ag1inr3 = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= sbkr[0] and group1.positions[:, dimension][i] <= sbkr[1])]
    ag2inr1 = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= sbkl[0] and group2.positions[:, dimension][i] < sbkl[1])]
    ag2inr2 = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= sbkm[0] and group2.positions[:, dimension][i] < sbkm[1])]
    ag2inr3 = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= sbkr[0] and group2.positions[:, dimension][i] <= sbkr[1])]
    if withsalt:
        group3 = uta.select_atoms("type " + ion3)
        ag3inr1 = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= sbkl[0] and group3.positions[:, dimension][i] < sbkl[1])]
        ag3inr2 = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= sbkm[0] and group3.positions[:, dimension][i] < sbkm[1])]
        ag3inr3 = [i for i in range(len(group3)) if (group3.positions[:, dimension][i] >= sbkr[0] and group3.positions[:, dimension][i] <= sbkr[1])]
        return ag1inr1, ag1inr2, ag1inr3, ag2inr1, ag2inr2, ag2inr3, ag3inr1, ag3inr2, ag3inr3
    else:
        return ag1inr1, ag1inr2, ag1inr3, ag2inr1, ag2inr2, ag2inr3


def rSubset(atomlist, r):
    # return list of all subsets of length r 
    # to deal with duplicate subsets use  
    # set(list(combinations(atomlist, r))) 
    return list(combinations(atomlist, r))

def solvation_shell(distpair):
    atoms_in_shell = []
    for i in range(len(distpair)):
        atoms_in_shell.append(distpair[i][1])
    return atoms_in_shell


def get_archs(raw_bridges, dp):
    bridges = 0
    for pair in raw_bridges:
        ion1, ion2 = pair
        t1 = int(math.floor(ion1 / dp) + 1)
        t2 = int(math.floor(ion2 / dp) + 1)
        if t1 != t2:
            bridges += 1
    return bridges


def get_bridges(top, trj, ion1, ion2, ion3, withsalt, dimension, concentration, cutoff1, cutoff2, dp, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
    ts = uta.trajectory[frame_id]
    cell = ts.dimensions
    if not withsalt:
        num_bridges = []
        for atom in range(len(group1)):
            distpair1 = distances.capped_distance(group1.positions[atom], group2.positions, cutoff1, box = cell)[0]
            atoms_in_shell = solvation_shell(distpair1)
            raw_bridges = rSubset(atoms_in_shell, 2)
            bridges = get_archs(raw_bridges, dp)
            num_bridges.append(bridges)
        return np.sum(num_bridges)
    else:
        type1, type2 = [], []
        type1_bridges, type2_bridges = [], []
        group3 = uta.select_atoms("type " + ion3)
        for atom in range(len(group1)):
            distpair1 = distances.capped_distance(group1.positions[atom], group2.positions, cutoff1, box = cell)[0]
            distpair2 = distances.capped_distance(group1.positions[atom], group3.positions, cutoff2, box = cell)[0]
            atoms_in_shell = solvation_shell(distpair1)
            raw_bridges = rSubset(atoms_in_shell, 2)
            bridges = get_archs(raw_bridges, dp)
            if len(distpair1) != 0 and len(distpair2) == 0:
                type1.append(atom)
                type1_bridges.append(bridges)
            elif len(distpair1) != 0 and len(distpair2) != 0:
                type2.append(atom)
                type2_bridges.append(bridges)
        return len(type1), len(type2), np.sum(type1_bridges), np.sum(type2_bridges)


def main():
    missions = [True]
    top = "../gentpr/cg_topol.tpr"
    trj = "../data_101.xtc"
    nt=20
    atom_group1 = "TF"
    atom_group2 = "IM"
    atom_group3 = "LI"
    withsalt = False
    dimension = 1  # 0: X; 1: Y; 2: Z
    concentration = "c000"
    cutoff1 = 7.7  # cutoff for anion - cation
    cutoff2 = 5.1  # cutoff for anion - lithium
    dp = 20  # conductive chain length
    uta = mda.Universe(top, trj)
    frame_ids = [ts.frame for ts in uta.trajectory]
    if missions[0]:
        do_analyse = functools.partial(get_bridges, top, trj, atom_group1, atom_group2, atom_group3, withsalt, dimension, concentration, cutoff1, cutoff2, dp)
        pool = mp.Pool(processes=nt)
        data = pool.imap(do_analyse, frame_ids)
        if withsalt:
            type1, type2 = [], []
            type1_bridges, type2_bridges = [], []
            for d1, d2, d3, d4 in data:
                type1.append(d1)
                type2.append(d2)
                type1_bridges.append(d3)
                type2_bridges.append(d4)
            bridgePrint1(concentration + '_bridge_pros', type1, type2, type1_bridges, type2_bridges)
        else:
            num_bridges = []
            for d1 in data:
                num_bridges.append(d1)
            bridgePrint2(concentration + '_bridge_pros', num_bridges)
            
if __name__ == "__main__":
    main()