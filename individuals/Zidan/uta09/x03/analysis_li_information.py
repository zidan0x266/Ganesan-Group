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


def get_li_infor(top, trj, ion1, ion2, ion3, DP, cutoff1, cutoff2, frame_id):
    uta = mda.Universe(top, trj)
    group1 = uta.select_atoms("type " + ion1)
    group2 = uta.select_atoms("type " + ion2)
    group3 = uta.select_atoms("type " + ion3)
    ts = uta.trajectory[frame_id]
    cell = ts.dimensions
    la, lc = [], []
    lch, lcoco = [], []
    for atom in range(len(group3)):
        nchains = []
        ncations = []
        distpair = distances.capped_distance(group3.positions[atom], group1.positions, cutoff2, box = cell)[0]
        la.append(len(distpair))
        a_id = [x[1] for x in distpair]
        for anion in a_id:
            distpair2 = distances.capped_distance(group1.positions[anion], group2.positions, cutoff1, box = cell)[0]
            c_id = [x[1] for x in distpair2]
            for cation in c_id:
                ncations.append(cation)
                chaindx = int(math.floor(cation / DP) + 1)
                nchains.append(chaindx)
        lc.append(len(set(ncations)))
        lcoco.append(len(ncations))
        lch.append(len(set(nchains)))
    return Counter(la), Counter(lc), Counter(lch), Counter(lcoco)


def main():
    top = "../cg_topol.tpr"
    trj = "data.xtc"
    nt=20
    atom_group1 = "TF"
    atom_group2 = "IM"
    atom_group3 = "LI"
    cutoff1 = 7.7  # cutoff for anion - cation
    cutoff2 = 5.1  # cutoff for anion - lithium
    dp = 20  # conductive chain length
    uta = mda.Universe(top, trj)
    frame_ids = [ts.frame for ts in uta.trajectory]
    freq_li_anion, freq_li_cation, freq_li_chain, freq_li_coco = [], [], [], []
    do_analyse = functools.partial(get_li_infor, top, trj, atom_group1, atom_group2, atom_group3, dp, cutoff1, cutoff2)
    pool = mp.Pool(processes=nt)
    data = pool.imap(do_analyse, frame_ids)
    for d1, d2, d3, d4 in data:
        freq_li_anion.append(d1)
        freq_li_cation.append(d2)
        freq_li_chain.append(d3)
        freq_li_coco.append(d4)
    li_anion_Num, li_anion_Pro = assoFreq(freq_li_anion)
    li_cation_Num, li_cation_Pro = assoFreq(freq_li_cation)
    li_chain_Num, li_chain_Pro = assoFreq(freq_li_chain)
    li_coco_Num, li_coco_Pro = assoFreq(freq_li_coco)
    assoPrint('li_anion', li_anion_Num, li_anion_Pro)
    assoPrint('li_cation', li_cation_Num, li_cation_Pro)
    assoPrint('li_chain', li_chain_Num, li_chain_Pro)
    assoPrint('li_coco', li_coco_Num, li_coco_Pro)
    print("Analysis of anion - cation associations complete!")


if __name__ == "__main__":
    main()
