#!/usr/bin/env

import mdtraj
import numpy as np
import scipy.integrate as integrate
from mdtraj.utils import ensure_type
from mdtraj.geometry.distance import compute_distances

#**************************************************************
#     this script computes the charge sum rule, namely
#       int_V[ 4*pi*r^2 (sum_a q_a*rho_a*g_ab(r))dr ] = -q_b
#       where 'a' and 'b' are ion types, e.g. either cation or anion
#       we use the 'un-normalized g(r)'s' to compute these quantities,
#       as normalizing in most RDF codes is not rigorously correct
#       as they don't subtract out origin ion in self-correlation,
#       which messes up the sum rules 
#*************************************************************

name_cat='name BB1'
name_an='name PFA'

nions=256
nions_2=128
bin_width=0.005
r_range=(0.0,2.2)
n_bins = int((r_range[1] - r_range[0]) / bin_width)

# load trajectory, topology
traj = mdtraj.load('mergextc/subdata.xtc', top = 'cg_traj.gro')
top=traj.topology

#****** first anion-anion
pairs=top.select_pairs(name_an,name_an)
distances = compute_distances(traj, pairs, periodic=True, opt=True)
# un-normalized g(r)
Gr_aa, edges = np.histogram(distances, range=r_range, bins=n_bins)
r = 0.5 * (edges[1:] + edges[:-1])

#****** now cation-anion
pairs=top.select_pairs(name_cat,name_an)
distances = compute_distances(traj, pairs, periodic=True, opt=True)
# un-normalized g(r)
Gr_ca, edges = np.histogram(distances, range=r_range, bins=n_bins)
r = 0.5 * (edges[1:] + edges[:-1])

#****** now cation-cation
pairs=top.select_pairs(name_cat,name_cat)
distances = compute_distances(traj, pairs, periodic=True, opt=True)
# un-normalized g(r)
Gr_cc, edges = np.histogram(distances, range=r_range, bins=n_bins)
r = 0.5 * (edges[1:] + edges[:-1])

# normalize by number of frames, and number of averaging centers per frame
# note for same-ion correlation, we only use unique pairs in compute_distances,
# so divide by nions_2
Gr_aa = Gr_aa / traj.n_frames / (nions_2)
Gr_cc = Gr_cc / traj.n_frames / (nions_2)
Gr_ca = Gr_ca / traj.n_frames / nions

gsum1 = np.subtract( Gr_ca , Gr_aa )
gsum2 = np.subtract( Gr_cc , Gr_ca )


for i in range(len(r)):
   # running integral
   q_x = np.sum(gsum1[:i])
   q_y = np.sum(gsum2[:i])
   print( r[i] , q_x , q_y )