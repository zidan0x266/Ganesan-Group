"""
Copyright (C) 2018-2019 Zidan Zhang <zhangzidan@gmail.com>
Jakub Krajniak <jkrajniak@gmail.com>

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

import math
from analysis import *


def extractHop(h5, dataset_name, counterN):
    """
    This function extracts hopping information.
    """
    ISAACS = rawLoad(h5, dataset_name)
    HOPES = []  # Chain association information
    hopes = []  # ion hopping information
    hop_ds = []
    hop_as = []
    for tf in range(len(ISAACS)):
        if tf <= len(ISAACS) - 2:
            association = []
            for i in range(len(ISAACS[tf])):
                anion, cation = ISAACS[tf][i]
                # 8 means everychain has 8 cations, should add a new argument later to define this by user.
                chain_id = int(math.floor(int(cation) / counterN) + 1)
                association.append((anion, chain_id))
            # fisrt remove the repeated pairs by set() then save the list format data
            HOPES.append(list(set(association)))
        if tf >= 1:
            hop_total = []
            hop_dissociation = []
            hop_association = []
            # Total hopping events between t and t - 1
            # DeltaH = list(set(ISAACS[tf - 1]) - set(ISAACS[tf])) + list(set(ISAACS[tf]) - set(ISAACS[tf - 1]))
            DeltaH = list(np.setdiff1d(ISAACS[tf - 1], ISAACS[tf])) + list(np.setdiff1d(ISAACS[tf], ISAACS[tf - 1]))
            # Dissociation events
            # DeltaDs = list(set(ISAACS[tf - 1]) - set(ISAACS[tf]))
            DeltaDs = list(np.setdiff1d(ISAACS[tf - 1], ISAACS[tf]))
            # Association events
            # DeltaAs = list(set(ISAACS[tf]) - set(ISAACS[tf - 1]))
            DeltaAs = list(np.setdiff1d(ISAACS[tf], ISAACS[tf - 1]))
            for i in range(len(DeltaH)):
                anion, cation = DeltaH[i]
                chain_id = int(math.floor(int(cation) / counterN) + 1)
                hop_total.append((anion, chain_id))
            hopes.append(hop_total)
            for i in range(len(DeltaDs)):
                anion, cation = DeltaDs[i]
                chain_id = int(math.floor(int(cation) / counterN) + 1)
                hop_dissociation.append((anion, chain_id))
            hop_ds.append(hop_dissociation)
            for i in range(len(DeltaAs)):
                anion, cation = DeltaAs[i]
                chain_id = int(math.floor(int(cation) / counterN) + 1)
                hop_association.append((anion, chain_id))
            hop_as.append(hop_association)
    return HOPES, hopes, hop_ds, hop_as


def analysis(h5, dataset_name, counterN):
    """
    This function analysis the hopping events.
    """
    HOPES, hopes, hop_ds, hop_as = extractHop(h5, dataset_name, counterN)
    intra_hop = []
    inter_hop = []
    intra_dissociation = []
    inter_dissociation = []
    intra_association = []
    inter_association = []
    for tf in range(len(HOPES)):
        inter_events = 0
        intra_events = 0
        for i in range(len(hopes[tf])):
            if hopes[tf][i] in HOPES[tf]:
                intra_events += 1
            else:
                inter_events += 1
        # intra_hop.append(float(intra_events) / (intra_events + inter_events))
        # inter_hop.append(float(inter_events) / (intra_events + inter_events))
        intra_hop.append(intra_events)
        inter_hop.append(inter_events)

        inter_ds = 0
        intra_ds = 0
        for i in range(len(hop_ds[tf])):
            if hop_ds[tf][i] in HOPES[tf]:
                intra_ds += 1
            else:
                inter_ds += 1
        # intra_dissociation.append(float(intra_ds) / (intra_ds + inter_ds))
        # inter_dissociation.append(float(inter_ds) / (intra_ds + inter_ds))
        intra_dissociation.append(intra_ds)
        inter_dissociation.append(inter_ds)

        inter_as = 0
        intra_as = 0
        for i in range(len(hop_as[tf])):
            if hop_as[tf][i] in HOPES[tf]:
                intra_as += 1
            else:
                inter_as += 1
        # intra_association.append(float(intra_as) / (intra_as + inter_as))
        # inter_association.append(float(inter_as) / (intra_as + inter_as))
        intra_association.append(intra_as)
        inter_association.append(inter_as)
    return [intra_hop, inter_hop], [intra_dissociation, inter_dissociation], [intra_association, inter_association]
