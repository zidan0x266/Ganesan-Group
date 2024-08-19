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


def ioninDPrint(filename, ioninD):
    with open('{}.xvg'.format(filename), 'w') as anaout:
        #print('# to be added', file=anaout)
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


def get_domains(curr_case):
    cases = ['LA10', 'LA15', 'LA20']
    if curr_case == cases[0]:
        domains = [14.18, 32.86, 66.80, 87.46]  # A30B10
    elif curr_case == cases[1]:
        domains = [22.06, 32.38, 83.64, 92.90]  # A45B15
    elif curr_case == cases[2]:
        domains = [25.42, 33.96, 93.06, 100.94]  # A60B20
    # define regions
    if1 = [domains[0], domains[1]]
    bk1 = [domains[1], domains[2]]
    if2 = [domains[2], domains[3]]
    to1 = [domains[0], domains[3]]
    return if1, bk1, if2, to1


def get_ions_inD(top, trj, ion1, ion2, dimension, curr_case, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
    ts = uta.trajectory[frame_id]
    if1, bk1, if2, to1 = get_domains(curr_case)
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
    data = [len(ag1into), len(ag1inbk), len(ag1inif), len(ag2into), len(ag2inbk), len(ag2inif)]
    return ts.frame, data


def get_ionlist_inD(top, trj, ion1, ion2, dimension, curr_case, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
    ts = uta.trajectory[frame_id]
    if1, bk1, if2, to1 = get_domains(curr_case)
    ag1inif1 = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= if1[0] and group1.positions[:, dimension][i] < if1[1])]
    ag1inif2 = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= if2[0] and group1.positions[:, dimension][i] <= if2[1])]
    ag1inbk = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= bk1[0] and group1.positions[:, dimension][i] < bk1[1])]
    ag1into = [i for i in range(len(group1)) if (group1.positions[:, dimension][i] >= to1[0] and group1.positions[:, dimension][i] <= to1[1])]
    #ag2inif1 = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= if1[0] and group2.positions[:, dimension][i] < if1[1])]
    #ag2inif2 = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= if2[0] and group2.positions[:, dimension][i] <= if2[1])]
    #ag2inbk = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= bk1[0] and group2.positions[:, dimension][i] < bk1[1])]
    #ag2into = [i for i in range(len(group2)) if (group2.positions[:, dimension][i] >= to1[0] and group2.positions[:, dimension][i] <= to1[1])]
    ag1inif = ag1inif1 + ag1inif2
    #ag2inif = ag2inif1 + ag2inif2
    #return ag1into, ag1inbk, ag1inif, ag2into, ag2inbk, ag2inif
    return ag1into, ag1inbk, ag1inif


def get_asso_domains(top, trj, ion1, ion2, dimension, curr_case, cutoff1, cutoff2, dp, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
    _, ag1inbk, ag1inif = get_ionlist_inD(top, trj, ion1, ion2, dimension, curr_case, frame_id)
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
        distpair = distances.capped_distance(group1.positions[atom], group2.positions, cutoff2, box = cell)[0]
        asso_pair_if.append(len(distpair))
        for i in range(len(distpair)):
            atom_pair = distpair[i][1]
            chain_ids.append(int(math.floor(atom_pair / dp) + 1))
        asso_chain_if.append(len(np.unique(chain_ids)))
    asso_pair_to = asso_pair_bk + asso_pair_if
    asso_chain_to = asso_chain_bk + asso_chain_if
    return Counter(asso_pair_to), Counter(asso_chain_to), Counter(asso_pair_bk), Counter(asso_chain_bk), Counter(asso_pair_if), Counter(asso_chain_if)


def main():
    missions_list = ["ions_in_domains", "ion_associations"]
    missions = [True, True]
    top = "../cg_topol.tpr"
    trj = "../cgtrj.xtc"
    nt=20
    atom_group1 = "TF"
    atom_group2 = "IM"
    dimension = 1  # 0: X; 1: Y; 2: Z
    curr_case = "LA10"
    dp = 10  # conductive chain length
    #cutoff1 = 8.0  # cutoff for anion - cation in Domain
    cutoff1 = 7.7  # cutoff for anion - cation in Bulk
    cutoff2 = 8.2  # cutoff for anion - cation in Interfacial region
    uta = mda.Universe(top, trj)
    frame_ids = [ts.frame for ts in uta.trajectory]
    if missions[0]:
        ioninD = []
        do_analyse = functools.partial(get_ions_inD, top, trj, atom_group1, atom_group2, dimension, curr_case)
        pool = mp.Pool(processes=nt)
        data = pool.imap(do_analyse, frame_ids)
        for _, ion_distribution in data:
            ioninD.append(ion_distribution)
        ioninDPrint(missions_list[0], ioninD)
        print("Analysis of ions in domains complete!")
    if missions[1]:
        freq_pair_to, freq_chain_to = [], []
        freq_pair_bk, freq_chain_bk = [], []
        freq_pair_if, freq_chain_if = [], []
        do_analyse = functools.partial(get_asso_domains, top, trj, atom_group1, atom_group2, dimension, curr_case, cutoff1, cutoff2, dp)
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


if __name__ == "__main__":
    main()