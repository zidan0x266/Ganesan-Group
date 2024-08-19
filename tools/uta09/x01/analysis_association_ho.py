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


def get_asso_12(top, trj, ion1, ion2, ion3, withsalt, cutoff1, cutoff2, dp, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
    if withsalt:
        group3 = uta.select_atoms("type " + ion3)
    ts = uta.trajectory[frame_id]
    cell = ts.dimensions
    asso_pair, asso_chain = [], []
    for atom in range(len(group1)):
        chain_ids = []
        distpair = distances.capped_distance(group1.positions[atom], group2.positions, cutoff1, box = cell)[0]
        asso_pair.append(len(distpair))
        for i in range(len(distpair)):
            atom_pair = distpair[i][1]
            chain_ids.append(int(math.floor(atom_pair / dp) + 1))
        asso_chain.append(len(np.unique(chain_ids)))
    return Counter(asso_pair), Counter(asso_chain)


def get_asso_13(top, trj, ion1, ion2, ion3, withsalt, cutoff1, cutoff2, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    if withsalt:
        group3 = uta.select_atoms("type " + ion3)
    ts = uta.trajectory[frame_id]
    cell = ts.dimensions
    asso_pair1, asso_pair2 = [], []
    for atom in range(len(group1)):
        distpair = distances.capped_distance(group1.positions[atom], group3.positions, cutoff2, box = cell)[0]
        asso_pair1.append(len(distpair))
    for atom in range(len(group3)):
        distpair = distances.capped_distance(group3.positions[atom], group1.positions, cutoff2, box = cell)[0]
        asso_pair2.append(len(distpair))
    return Counter(asso_pair1), Counter(asso_pair2)


def main():
    missions = [True, False]
    top = "../gentpr/cg_topol.tpr"
    trj = "../data_101.xtc"
    nt=20
    atom_group1 = "TF"
    atom_group2 = "IM"
    atom_group3 = "LI"
    withsalt = True
    if withsalt:
        missions[1] = True
    cutoff1 = 7.7  # cutoff for anion - cation
    cutoff2 = 5.1  # cutoff for anion - lithium
    dp = 20  # conductive chain length
    uta = mda.Universe(top, trj)
    frame_ids = [ts.frame for ts in uta.trajectory]
    if missions[0]:
        freq_pair, freq_chain = [], []
        do_analyse = functools.partial(get_asso_12, top, trj, atom_group1, atom_group2, atom_group3, withsalt, cutoff1, cutoff2, dp)
        pool = mp.Pool(processes=nt)
        data = pool.imap(do_analyse, frame_ids)
        for d1, d2 in data:
            freq_pair.append(d1)
            freq_chain.append(d2)
        pair_Num, pair_Pro = assoFreq(freq_pair)
        chain_Num, chain_Pro = assoFreq(freq_chain)
        assoPrint('anion_cation' + '_Pair', pair_Num, pair_Pro)
        assoPrint('anion_cation' + '_Chain', chain_Num, chain_Pro)
        print("Analysis of anion - cation associations complete!")
    if missions[1]:
        freq_pair_13, freq_pair_31 = [], []
        do_analyse = functools.partial(get_asso_13, top, trj, atom_group1, atom_group2, atom_group3, withsalt, cutoff1, cutoff2)
        pool = mp.Pool(processes=nt)
        data = pool.imap(do_analyse, frame_ids)
        for d1, d2 in data:
            freq_pair_13.append(d1)
            freq_pair_31.append(d2)
        pair_Num_13, pair_Pro_13 = assoFreq(freq_pair_13)
        pair_Num_31, pair_Pro_31 = assoFreq(freq_pair_31)
        assoPrint('anion_lithium' + '_Pair', pair_Num_13, pair_Pro_13)
        assoPrint('lithium_anion' + '_Pair', pair_Num_31, pair_Pro_31)
        print("Analysis of anion - lithium associations complete!")


if __name__ == "__main__":
    main()