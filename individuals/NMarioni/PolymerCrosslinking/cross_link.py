# This script is designed to crosslink MPD and TMC monomers to create a crosslinked polyamide membrane
# This script can be adapted to crosslinked different membranes

import builder_io # Provided in the mdp folder
import copy
import numpy as np
import h5py
import os
import MDAnalysis as mda
import MDAnalysis.lib.distances as distances
import networkx as nx
import time

from sys import argv
script, gro_file, out_gro, top_file, out_top = argv

# This is a topology file after crosslinking the two monomers, here MPD and TMC
Charge_TOP = builder_io.GroTopFile('mdp/MPDTMC.top'); Charge_TOP.read()

type_arr = []; charge_arr = []
for i in Charge_TOP.atoms:
    type_arr.append(Charge_TOP.atoms[i].atom_type); charge_arr.append(Charge_TOP.atoms[i].charge)
type_arr, bool_arr = np.unique(np.array(type_arr), return_index=True); charge_arr = np.array(charge_arr)[bool_arr]





# Molecule residue names; molecule crosslinker names, molecule crosslinker names for MDAnalysis
res_a1 = 'MPD'; names_a1 = np.array(['NA1','NA2']); uta_a1 = 'name NA1 or name NA2'
res_a2 = 'TMC'; names_a2 = np.array(['CC3','CC4','CC5']); uta_a2 = 'name CC3 or name CC4 or name CC5'
# Defines alterations to be made to crosslinker atom + associated atoms
# Defines a tuple per crosslinker atom
# tuple[0] = define crosslinker atom to be updated
#            crosslinker atom name, new atom type, new atom name, e.g., 'NA1 opls_238 NA' updates crosslinker NA1 to have atom type = opls_238, atom name = NA
# tuple[1:x] = define atoms to be updated -> starting with atoms that appear earliest in the residue
#              atom name, instance of that name, new atom type, new atom name, e.g., 'HA1 2 opls_241 HA' updates the second instance of HA1, starting from the end of the residue, to have atom type = opls_241, atom name = HA
# tuple[x:] = define atoms to be deleted, relative to the END of the residue -> starting with atoms that appear latest in the residue
#             del, atom name, instance of that name, e.g., 'del HA1 1' deletes the first instance of HA1, starting from the end of the residue
alter_a1 = [('NA1 opls_238 NA','CB1 1 opls_230 CB','HA1 2 opls_241 HA','del HA1 1'),('NA2 opls_238 NA','CB2 1 opls_230 CB','HA2 2 opls_241 HA','del HA2 1')]
alter_a2 = [('CC3 opls_234 CC','CB3 1 opls_145 CB','OC3 2 opls_236 OC','del HC3 1','del OC3 1'),('CC4 opls_234 CC','CB4 1 opls_145 CB','OC4 2 opls_236 OC','del HC4 1','del OC4 1'),('CC5 opls_234 CC','CB5 1 opls_145 CB','OC5 2 opls_236 OC','del HC5 1','del OC5 1')]
# Define improper dihedrals for crosslinker a1 and a2, if one is formed upon crosslinking, e.g., improper_a1 = ['1','improper_Z_N_X_Y'] for improper_Z_N_X_Y, function 1 
# If imrpoper dihedral is not formed, first position should be '-1'
improper_a1 = ['1','improper_Z_N_X_Y']; improper_a2 = ['1','improper_O_C_X_Y']

## If files exist, delete them
#if os.path.isfile(out_top):
#    os.remove(out_top)
#if os.path.isfile(out_gro):
#    os.remove(out_gro)



### MDAnalysis code to find bonds ###
uta = mda.Universe(gro_file)

# Extract crosslinker atom positions and cell size
ts = uta.trajectory[0]
a1 = uta.select_atoms(uta_a1).positions; n_a1 = len(a1)
a2 = uta.select_atoms(uta_a2).positions; n_a2 = len(a2)
cell = ts.dimensions; cell[2] *= 10 # Make z-axis longer to preventing crosslinking across the z axis

# Calculate distances between crosslinkers
pair, dist = distances.capped_distance(a1, a2, 4.5, box=cell); del a1; del a2

print("Current Crosslink Density")
print(1 - (n_a1/1380.0))

Finalize = False
# Define the desired degree of crosslinking, where (1 - (n_a1/1380.0)) = 0.90 defines that point at 90% of possible crosslinks are formed
if (1.0 - (n_a1/1380.0)) >= 0.90:
    print("Desired Crosslink Density Achieved")
    print(1 - (n_a1/1380.0))
    Finalize = True

    # The code below is to track the Crosslinking Process paired with the code described above
    count = 50
    with open('xlink_tracker.txt','w') as f:
        f.write(str(count))
elif len(pair) == 0:
    print("No Crosslinked Formed")



    # The code below is to track the Crosslinking Process if it is desired to stop the process after not forming any new crosslinks {count} times in a row
    if os.path.isfile('xlink_tracker.txt'):
        with open('xlink_tracker.txt','r') as f:
            f = f.readlines()
        count = int(f[0]) + 1
    else:
        count = 1
    
    with open('xlink_tracker.txt','w') as f:
        f.write(str(count))
    
    if count == 50:
        Finalize = True



    if Finalize == False:
        exit()
elif os.path.isfile('xlink_tracker.txt'):
    os.remove('xlink_tracker.txt')

if Finalize:
    # This if statement deletes un-crosslinked monomer. This can be monomer that forms no crosslinks, or polymer chains that are not crosslinked to the continuous membrane network
    # If you do not want any monomer deleted, comment out everything inside this if statement except for exit()



    # This section finds all monomer that formed no crosslinks for removal from the system
    track_a1 = []; tracker = 0; counter = 0
    for i in uta.select_atoms(uta_a1):
        i = i.resid
        if i != tracker:
            tracker = i; counter = 1
        else:
            counter += 1
        
        if counter == len(names_a1):
            track_a1.append(i)
        elif counter > len(names_a1):
            print("Error in Processing Number of Non-Xlinked Molecules")
            exit()
    
    track_a2 = []; tracker = 0; counter = 0
    for i in uta.select_atoms(uta_a2):
        i = i.resid
        if i != tracker:
            tracker = i; counter = 1
        else:
            counter += 1
        
        if counter == len(names_a2):
            track_a2.append(i)
        elif counter > len(names_a2):
            print("Error in Processing Number of Non-Xlinked Molecules")
            exit()
    
    print('Type 1 Molecules with no xlinks: ')
    print(track_a1)
    print('Type 2 Molecules with no xlinks: ')
    print(track_a2)
    print('Removing Molecules')



    # Read in GRO and TOP files to be edited
    GRO = builder_io.GROFile(gro_file); GRO.read()
    TOP = builder_io.GroTopFile(top_file); TOP.read()



    # This section find all atoms that are continuously bonded to each other. All atoms that are not a part of the largest network (the continuous membrane) are deleted.
    # This section makes the above section for un-crosslinked monomers redundant
    track_a1 = []; track_a2 = []

    G = nx.Graph(); G.add_nodes_from(range(1,len(TOP.atoms)))
    for a in TOP.bonds_def:
        for b in TOP.bonds_def[a]:
            G.add_edge(a,b)
    
    for i, z in enumerate(sorted(nx.connected_components(G), key = len, reverse=True)):
        z = np.array(list(z))

        # The first group is the main crosslinked membrane -> DO NOT DELETE ATOMS
        if i == 0:
            print("Number of Atoms in the continuous membrane: " + str(len(z)))
            continue

        print("Number of Atoms in monomer/polymer cluster: " + str(len(z)))
        for z_j in z:
            if TOP.atoms[z_j].chain_name == res_a1:
                if TOP.atoms[z_j].chain_idx in track_a1:
                    continue
                track_a1.append(TOP.atoms[z_j].chain_idx)
            elif TOP.atoms[z_j].chain_name == res_a2:
                if TOP.atoms[z_j].chain_idx in track_a2:
                    continue
                track_a2.append(TOP.atoms[z_j].chain_idx)
    
    print('Type 1 Molecules not crosslinked to the continuous membrane: ')
    print(track_a1)
    print('Type 2 Molecules not crosslinked to the continuous membrane: ')
    print(track_a2)
    print('Removing Molecules')



    ### Update GRO and TOP atoms section -> remove atoms ###
    # Remove atoms
    n = len(GRO.atoms); deleted = []
    res_save = 0; res_index = 0
    for i in range(n, 0, -1):

        if GRO.atoms[i].chain_idx in track_a1 or GRO.atoms[i].chain_idx in track_a2:
            deleted.append(i)
            del GRO.atoms[i]; del TOP.atoms[i]
    deleted = np.array(deleted); min_del = min(deleted)
    # Renumber atoms to account for deletions
    GRO.renumber(); TOP.renumber()

    ### Renumber Atoms, Bonds, Angles, Dihedrals, Improper Dihedrals, Pairs ###
    # Bonds
    for b in copy.copy(TOP.bonds):
        if b[0] in deleted or b[1] in deleted:
            del TOP.bonds[b]
        elif b[1] > min_del or b[0] > min_del:
            bt = copy.copy(TOP.bonds[b])
            del TOP.bonds[b]
            TOP.bonds[(b[0] - len(deleted[deleted < b[0]]), b[1] - len(deleted[deleted < b[1]]))] = bt
    del b; del bt
    # Angles
    for a in copy.copy(TOP.angles):
        if a[0] in deleted or a[1] in deleted or a[2] in deleted:
            del TOP.angles[a]
        elif a[2] > min_del or a[1] > min_del or a[0] > min_del:
            at = copy.copy(TOP.angles[a])
            del TOP.angles[a]
            TOP.angles[(a[0] - len(deleted[deleted < a[0]]), a[1] - len(deleted[deleted < a[1]]), a[2] - len(deleted[deleted < a[2]]))] = at
    del a; del at
    # Dihderals
    for d in copy.copy(TOP.dihedrals):
        if d[0] in deleted or d[1] in deleted or d[2] in deleted or d[3] in deleted:
            del TOP.dihedrals[d]
        elif d[3] > min_del or d[2] > min_del or d[1] > min_del or d[0] > min_del:
            dt = copy.copy(TOP.dihedrals[d])
            del TOP.dihedrals[d]
            TOP.dihedrals[(d[0] - len(deleted[deleted < d[0]]), d[1] - len(deleted[deleted < d[1]]), d[2] - len(deleted[deleted < d[2]]), d[3] - len(deleted[deleted < d[3]]))] = dt
    del d; del dt
    # Improper Dihedrals
    for i in copy.copy(TOP.improper_dihedrals):
        if i[0] in deleted or i[1] in deleted or i[2] in deleted or i[3] in deleted:
            del TOP.improper_dihedrals[i]
        elif i[3] > min_del or i[2] > min_del or i[1] > min_del or i[0] > min_del:
            it = copy.copy(TOP.improper_dihedrals[i])
            del TOP.improper_dihedrals[i]
            TOP.improper_dihedrals[(i[0] - len(deleted[deleted < i[0]]), i[1] - len(deleted[deleted < i[1]]), i[2] - len(deleted[deleted < i[2]]), i[3] - len(deleted[deleted < i[3]]))] = it
    del i; del it
    # Pairs
    for p in copy.copy(TOP.pairs):
        if p[0] in deleted or p[1] in deleted:
            del TOP.pairs[p]
        elif p[1] > min_del or p[0] > min_del:
            pt = copy.copy(TOP.pairs[p])
            del TOP.pairs[p]
            TOP.pairs[(p[0] - len(deleted[deleted < p[0]]), p[1] - len(deleted[deleted < p[1]]))] = pt
    del p; del pt

    GRO.write(out_gro); TOP.write(out_top)



    exit()

# Two crosslinkers from a1 may try to bond with the same crosslinker of a2, or vice versa. Only bond the closest pair of crosslinkers
print("Number of Pairs (Including Multi-Pairs): " + str(len(pair)))
pair_arr = np.zeros((n_a1, n_a2))
for i in range(len(pair)):
    pair_arr[pair[i,0], pair[i,1]] = dist[i] # 2D array of distance between crosslinker_a1 i and crosslinker_a2 j
del pair; del dist

# If 2 crosslinker_a2 bonded to crosslinker_a1, only keep closest pair
pair_temp = np.zeros(n_a1, dtype=int) - 1
for i in range(n_a1):
    pair_i = pair_arr[i,:]
    if len(pair_i[pair_i != 0]) > 0:
        pair_temp[i] = np.where(pair_i == min(pair_i[pair_i != 0]))[0][0]

# If 2 crosslinker_a1 bonded to crosslinker_a2, only keep closest pair
pair_a1 = np.zeros(n_a1, dtype=int) - 1
for j in range(n_a2):
    pair_j = np.where(pair_temp == j)[0]
    if len(pair_j) > 0:
        pair_a1[np.where(pair_arr == min(pair_arr[pair_j, j]))[0][0]] = j
print("Number of Pairs (Excluding Multi-Pairs): " + str(len(pair_a1[pair_a1 != -1])))

# Limits molecule 1 to only crosslinking a single site per interation
for i in range(n_a1):
    if i%len(names_a1) != 0:
        continue
    pair_temp = np.array(pair_a1[i:i+len(names_a1)])

    if len(pair_temp[pair_temp != -1]) > 1:
        index_save = -1; dist_save = 0
        for j,k in enumerate(pair_temp):
            if k == -1:
                continue

            if index_save == -1:
                index_save = i+j; dist_save = pair_arr[i+j,k]
            elif pair_arr[i+j,k] == 0:
                print('Broken')
                exit()
            elif pair_arr[i+j,k] < dist_save:
                pair_a1[index_save] = -1
                index_save = i+j; dist_save = pair_arr[i+j,k]
            else:
                pair_a1[i+j] = -1
print("Number of Crosslinks (1-per a1 Molecule): " + str(len(pair_a1[pair_a1 != -1])))

# Limits molecule 2 to only crosslinking a single site per interation
for i in range(n_a2):
    if i%len(names_a2) != 0:
        continue

    pair_temp  = np.zeros(len(names_a2), dtype=int) - 1
    for j in range(len(names_a2)):
        temp = np.where(pair_a1 == i+j)[0]
        if len(temp) != 0:
            pair_temp[j] = temp[0]
    pair_temp = np.array(pair_temp)

    if len(pair_temp[pair_temp != -1]) > 1:
        index_save = -1; dist_save = 0
        for j,k in enumerate(pair_temp):
            if k == -1:
                continue
            
            if index_save == -1:
                index_save = k; dist_save = pair_arr[k,i+j]
            elif pair_arr[k,i+j] == 0:
                print('Broken')
                exit()
            elif pair_arr[k,i+j] < dist_save:
                pair_a1[index_save] = -1
                index_save = k; dist_save = pair_arr[k,i+j]
            else:
                pair_a1[k] = -1
print("Number of Crosslinks (1-per a2 Molecule): " + str(len(pair_a1[pair_a1 != -1])))



print("Number of Crosslinks (Final): " + str(len(pair_a1[pair_a1 != -1])))
del pair_arr; del pair_temp



### Update GRO and TOP atoms section -> remove atoms, update atoms types, charges ###
# Read in GRO and TOP files to be edited
GRO = builder_io.GROFile(gro_file); GRO.read()
TOP = builder_io.GroTopFile(top_file); TOP.read()

# Remove atoms, updates atom types, charges
n = len(GRO.atoms); deleted = []
count_a1 = n_a1; count_a2 = n_a2; res_save = 0; res_index = 0
for i in range(n, 0, -1):
    if i%1000 == 0:
        print('Remaining Crosslinkers: ' + str(count_a1 + count_a2))

    # Track index of last atom in current residue
    if GRO.atoms[i].chain_idx != res_save:
        res_save = GRO.atoms[i].chain_idx
        res_index = i
    
    # If crosslink formed, retrieve correct alteration tuple
    if GRO.atoms[i].name in names_a1 and pair_a1[count_a1-1] != -1:
        alter = alter_a1[np.where(GRO.atoms[i].name == names_a1)[0][0]]
        count_a1 -= 1
    elif GRO.atoms[i].name in names_a2 and count_a2-1 in pair_a1:
        alter = alter_a2[np.where(GRO.atoms[i].name == names_a2)[0][0]]
        count_a2 -= 1
    else:
        if GRO.atoms[i].name in names_a1:
            count_a1 -= 1
        elif GRO.atoms[i].name in names_a2:
            count_a2 -= 1
        continue

    # Perform alterations
    check_del = 0
    for k, alt in enumerate(alter):
        alt = alt.split(' ')

        # Edit Crosslinker Atom
        if k == 0:
            # Do not rename GRO.atoms[i].name of crosslinker atom yet
            TOP.atoms[i].atom_type = alt[1]; TOP.atoms[i].name = alt[2]
            TOP.atoms[i].charge = charge_arr[np.where(TOP.atoms[i].atom_type == type_arr)[0][0]]
        # Delete atoms
        elif alt[0] == 'del':
            if check_del == 0:
                n_save = res_index; check_del = 1
            count_temp = 0
            for n_temp in range(n_save, 0, -1):
                if n_temp in deleted:
                    continue
                if GRO.atoms[n_temp].name == alt[1]:
                    count_temp += 1
                if count_temp == int(alt[2]):
                    deleted.append(n_temp)
                    del GRO.atoms[n_temp]; del TOP.atoms[n_temp]
                    n_save = n_temp-1; break
        # Edit atoms
        else:
            count_temp = 0
            for n_temp in range(res_index, 0, -1):
                if n_temp in deleted:
                    continue
                if GRO.atoms[n_temp].name == alt[0]:
                    count_temp += 1
                if count_temp == int(alt[1]):
                    GRO.atoms.update({n_temp: GRO.atoms[n_temp]._replace(name = alt[3])})
                    TOP.atoms[n_temp].atom_type = alt[2]; TOP.atoms[n_temp].name = alt[3]
                    TOP.atoms[n_temp].charge = charge_arr[np.where(TOP.atoms[n_temp].atom_type == type_arr)[0][0]]
                    break
deleted = np.array(deleted); min_del = min(deleted)
# Renumber atoms to account for deletions
GRO.renumber(); TOP.renumber()

#copy.copy(GRO).write('removed_atoms.gro')
#copy.copy(TOP).write('removed_atoms.top')



### Renumber Atoms, Bonds, Angles, Dihedrals, Improper Dihedrals, Pairs ###
# Bonds
for b in copy.copy(TOP.bonds):
    if b[0] in deleted or b[1] in deleted:
        del TOP.bonds[b]
    elif b[1] > min_del or b[0] > min_del:
        bt = copy.copy(TOP.bonds[b])
        del TOP.bonds[b]
        TOP.bonds[(b[0] - len(deleted[deleted < b[0]]), b[1] - len(deleted[deleted < b[1]]))] = bt
del b; del bt
# Angles
for a in copy.copy(TOP.angles):
    if a[0] in deleted or a[1] in deleted or a[2] in deleted:
        del TOP.angles[a]
    elif a[2] > min_del or a[1] > min_del or a[0] > min_del:
        at = copy.copy(TOP.angles[a])
        del TOP.angles[a]
        TOP.angles[(a[0] - len(deleted[deleted < a[0]]), a[1] - len(deleted[deleted < a[1]]), a[2] - len(deleted[deleted < a[2]]))] = at
del a; del at
# Dihderals
for d in copy.copy(TOP.dihedrals):
    if d[0] in deleted or d[1] in deleted or d[2] in deleted or d[3] in deleted:
        del TOP.dihedrals[d]
    elif d[3] > min_del or d[2] > min_del or d[1] > min_del or d[0] > min_del:
        dt = copy.copy(TOP.dihedrals[d])
        del TOP.dihedrals[d]
        TOP.dihedrals[(d[0] - len(deleted[deleted < d[0]]), d[1] - len(deleted[deleted < d[1]]), d[2] - len(deleted[deleted < d[2]]), d[3] - len(deleted[deleted < d[3]]))] = dt
del d; del dt
# Improper Dihedrals
for i in copy.copy(TOP.improper_dihedrals):
    if i[0] in deleted or i[1] in deleted or i[2] in deleted or i[3] in deleted:
        del TOP.improper_dihedrals[i]
    elif i[3] > min_del or i[2] > min_del or i[1] > min_del or i[0] > min_del:
        it = copy.copy(TOP.improper_dihedrals[i])
        del TOP.improper_dihedrals[i]
        TOP.improper_dihedrals[(i[0] - len(deleted[deleted < i[0]]), i[1] - len(deleted[deleted < i[1]]), i[2] - len(deleted[deleted < i[2]]), i[3] - len(deleted[deleted < i[3]]))] = it
del i; del it
# Pairs
for p in copy.copy(TOP.pairs):
    if p[0] in deleted or p[1] in deleted:
        del TOP.pairs[p]
    elif p[1] > min_del or p[0] > min_del:
        pt = copy.copy(TOP.pairs[p])
        del TOP.pairs[p]
        TOP.pairs[(p[0] - len(deleted[deleted < p[0]]), p[1] - len(deleted[deleted < p[1]]))] = pt
del p; del pt

#copy.copy(TOP).write('renumbered.top')



### Re-read in TOP with renumbered bonds_def ###
TOP.write(out_top); TOP = builder_io.GroTopFile(out_top); TOP.read(); os.remove(out_top)



### Add new Bonds, Angles, Dihedrals, Improper Dihedrals, Pairs for crosslinks ###
# Find new crosslinker atom indices
n = len(GRO.atoms)
loc_a1 = np.zeros(n_a1, dtype=int); loc_a2 = np.zeros(n_a1, dtype=int)
count_a1 = 0; count_a2 = 0
for i in range(1, n + 1):
    if GRO.atoms[i].name in names_a1:
        loc_a1[count_a1] = i; count_a1 += 1
    elif GRO.atoms[i].name in names_a2:
        loc_a2[count_a2] = i; count_a2 += 1

# Create new bonds, angles, dihedrals, impropers, and pairs
for i,j in enumerate(pair_a1):
    if j == -1:
        continue

    # Now we can re-name the bonded crosslinkers in the GRO file
    GRO.atoms.update({loc_a1[i]: GRO.atoms[loc_a1[i]]._replace(name = TOP.atoms[loc_a1[i]].name)})
    GRO.atoms.update({loc_a2[j]: GRO.atoms[loc_a2[j]]._replace(name = TOP.atoms[loc_a2[j]].name)})

    G = nx.Graph()

    # Crosslinker a1 bonding list
    l_a1 = []; l_a1.append(loc_a1[i])
    for k in copy.copy(TOP.bonds_def[l_a1[0]]):
        l_a1.append(k); G.add_edge(l_a1[0],k)
        for l in copy.copy(TOP.bonds_def[k]):
            l_a1.append(l); G.add_edge(k,l)
    l_a1 = np.unique(l_a1)
    # Crosslinker a2 bonding list
    l_a2 = []; l_a2.append(loc_a2[j])
    for k in copy.copy(TOP.bonds_def[l_a2[0]]):
        l_a2.append(k); G.add_edge(l_a2[0],k)
        for l in copy.copy(TOP.bonds_def[k]):
            l_a2.append(l); G.add_edge(k,l)
    l_a2 = np.unique(l_a2)
    # Crosslink a1 and a2
    G.add_edge(loc_a1[i],loc_a2[j])
    
    # New bonds, angles, dihedrals, pairs
    # Define all bonding combinations (1,2,3,4,5,6), where 3,4 are the crosslinker atoms. This way, we cover all new interactions between the newly bonded molecules:
    # Bonds (3,4); Angles (2,3,4), (3,4,5); Dihedrals (1,2,3,4), (2,3,4,5), (3,4,5,6); Pairs (1,4), (2,5), (3,6), and any necessary improper dihedrals.
    # Note, for exmaple, that Angles (1,2,3) are not needed, as they are contianed within an exisiting molecule, etc
    for k in l_a1:
        for path in nx.all_simple_paths(G, source=k, target=l_a2, cutoff=3):
            if len(path) == 2:
                TOP.bonds[tuple(path)] = copy.copy(TOP.bonds[tuple((1,2))])
            elif len(path) == 3:
                TOP.angles[tuple(path)] = copy.copy(TOP.angles[tuple((1,2,3))])
            elif len(path) == 4:
                TOP.dihedrals[tuple(path)] = copy.copy(TOP.dihedrals[tuple((1,2,3,4))])
                TOP.pairs[tuple((path[0],path[3]))] = copy.copy(TOP.pairs[tuple((1,4))])
            else:
                print("Error in NetworkX Pathing")
                exit()
    # New improper dihedrals -> may need custom code for any extra dihedrals that do not involved crosslinker a1 or a2 as the central atom
    if improper_a1[0] != '-1':
        i_a1 = list(copy.copy(TOP.bonds_def[loc_a1[i]]))
        TOP.improper_dihedrals[(i_a1[0],i_a1[1],loc_a1[i],loc_a2[j])] = improper_a1
    if improper_a2[0] != '-1':
        i_a2 = list(copy.copy(TOP.bonds_def[loc_a2[j]]))
        TOP.improper_dihedrals[(i_a2[0],i_a2[1],loc_a2[j],loc_a1[i])] = improper_a2

#copy.copy(GRO).write('renamed.gro')
#copy.copy(TOP).write('crosslinked.top')



### Write New GRO and TOP ###
GRO.write(out_gro)
TOP.write(out_top)
