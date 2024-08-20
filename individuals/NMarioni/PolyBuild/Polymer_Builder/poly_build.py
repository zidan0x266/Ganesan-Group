import numpy as np
import h5py
import random
import os

from sys import argv

script, poly_names, HT_name, num_units, charge_scale, folder, sample, poly_type = argv
poly_names = np.array(poly_names.split(), dtype=np.str_); num_units = np.array(num_units.split(), dtype=int); charge_scale = float(charge_scale); poly_type = int(poly_type)
#poly_names = array of polymer names to be made into statistical copolymer
#num_units = array of the number of units per poly_names
#poly_type = type of polymer; 0 -> Random Statisical Polymer, 1 -> Block Copolymer
#e.g. poly_names = [PEG, LGA], num_units = [10,10] makes a 20mer of 50% PEG 50% LGA
#
#folder = folder to contain the sample_x folders; sample = sample number (for creating easy folders for making systems)

b_Tot = sum(num_units)
names = ['h', 'm', 't', 'hm', 'mt', 'hmt']

# If rotate_check == 1, then the linear polymer.gro file will have randomly rotated units
# If rotate_check == 0, then there will be no rotation
#       NOTE: Make sure that rotate works properly before using it
rotate_check = 1



def writeFile(name, array):
    f = open(name, 'w')
    for i in array:
        f.write(i)
    f.close()



def PolySeq_Rand():
# Determine the sequence of polymer units randomly based on the mol % per unit

    # Create array of poly_names for head and tail units
    head_names = []; tail_names = []
    for n in poly_names:
        head_names.append(n+'H')
        tail_names.append(n+'T')

    block = num_units*0                                     # array to track how many of each poly unit is added to the sequence
    Copolymer = np.empty(b_Tot, dtype = (np.unicode_, 16))  # array of poly unit names
    Copoly = np.empty(b_Tot, dtype = int)                   # array of poly unit indices (where indices is position in poly_names)

    if HT_name == 'None' or HT_name == '':
        count = 0
    else:
        count = 1; block[np.where(poly_names == HT_name)[0][0]] += 2
        Copolymer[0] = HT_name+'H'; Copolymer[-1] = HT_name+'T'
        Copoly[0] = np.where(poly_names == HT_name)[0][0]; Copoly[-1] = np.where(poly_names == HT_name)[0][0]

    while sum(block) < b_Tot:
        probs = (num_units - block) / (b_Tot - sum(block))  #array of probabilities that the next unit is poly_names[i]
        #print(probs)

        p = random.random() # random number between 0 and 1
        for i, p_i in enumerate(np.cumsum(probs)):
            if p < p_i:
                if count == 0:
                    Copolymer[count] = head_names[i]

                elif count == b_Tot-1:
                    Copolymer[count] = tail_names[i]

                else:
                    Copolymer[count] = poly_names[i]

                Copoly[count] = i
                block[i] += 1

                break
            
        count += 1

    print(Copolymer)

    return Copoly



def PolySeq_Block():
# Determine the sequence of polymer units, Block Copolymer

    # Create array of poly_names for head and tail units
    head_names = []; tail_names = []
    for n in poly_names:
        head_names.append(n+'H')
        tail_names.append(n+'T')

    Copolymer = np.empty(b_Tot, dtype = (np.unicode_, 16))  # array of poly unit names
    Copoly = np.empty(b_Tot, dtype = int)                   # array of poly unit indices (where indices is position in poly_names)

    count = 0
    for i, name in enumerate(poly_names):
        for j in range(num_units[i]):
            if count == 0:
                Copolymer[count] = head_names[i]
            elif count == b_Tot-1:
                Copolymer[count] = tail_names[i]
            else:
                Copolymer[count] = poly_names[i]
            
            Copoly[count] = i
            count += 1

    print(Copolymer)

    return Copoly



def PolyRead():
# Read in the h5py files containing the chopped up gro and topol files for each polymer created with split_top.py

    mon_len = []; mon_nums = []; head = []; tail = []
    gro_ar = []; atoms_ar = []; bonds_ar = []; angles_ar = []; dihedrals_ar = []; impropers_ar = []; pairs_ar = []
    for poly in poly_names:

        # split_top.py is within MKTOP, which contains a folder for each polymer unit with the .hdf5 file
        with h5py.File('../MKTOP/'+poly+'/top.hdf5','r') as f:
            mon_len.append(f['mon_len'][:])
            mon_nums.append(f['mon_nums'][:])
            head.append(f['head'].asstr()[:])
            tail.append(f['tail'].asstr()[:])

            gro = []; atoms = []; bonds = []; angles = []; dihedrals = []; impropers = []; pairs = []
            for i, n in enumerate(names):
                if i< 3:
                    gro.append(f[n+"_gro"].asstr()[:])
                    atoms.append(f[n+"_atoms"].asstr()[:])
                bonds.append(f[n+"_bonds"].asstr()[:])
                angles.append(f[n+"_angles"].asstr()[:])
                dihedrals.append(f[n+"_dihedrals"].asstr()[:])
                impropers.append(f[n+"_improper"].asstr()[:])
                pairs.append(f[n+"_pairs"].asstr()[:])
            gro_ar.append(gro)
            atoms_ar.append(atoms)
            bonds_ar.append(bonds)
            angles_ar.append(angles)
            dihedrals_ar.append(dihedrals)
            impropers_ar.append(impropers)
            pairs_ar.append(pairs)

    return mon_len, mon_nums, head, gro_ar, atoms_ar, bonds_ar, angles_ar, dihedrals_ar, impropers_ar, pairs_ar, tail



def editGro(gro, atom_num, delta, count):
# Edits lines in the .gro file to update the atom number, residue number, and position (to create a linear chain in the z-direction)
#
# Inputs: gro = array of .gro lines for the polymer unit; atom_num = number of atoms to be added so the numbers are sequential across the whole file;
#         delta = value to be added to the z axis; count = number of residues to be added so the numbers are sequential across the whole file

    # Randomly rotates the monomer unit around the z axis
    rand = random.random()
    if rotate_check == 0:
        rotate = [ 1, 1]
    elif rand < 0.25:
        rotate = [ 1, 1]
    elif rand < 0.50:
        rotate = [-1, 1]
    elif rand < 0.75:
        rotate = [-1,-1]
    else:
        rotate = [ 1,-1]

    edGro  = []
    for i, line in enumerate(gro):
        res_nr = int(line[:5]); nr = int(line.split()[2]); pos = np.array(line.split()[3:], dtype = float)
        edGro.append('{:5d}'.format(res_nr+count) + line[5:15] + '{:5d} {: 7.3f} {: 7.3f} {: 7.3f}'.format(nr+atom_num, pos[0]*rotate[0], pos[1]*rotate[1], pos[2] + (delta)) + '\n')
    
    return edGro

def MakeGro(Copoly, gro, mon_len):
# Assembles the .gro file for the polymer
#
# Inputs: Copoly = sequence of poly unit indices, where the indices is related to poly_names; gro = array of .gro chunks for [poly_names,head mid tail];
#         mon_len = array of lengths of each unit monomer from poly_names

    # Determine the length of each monomer unit along the z-axis
    delta_ar = []
    for gro_i in gro:
        min_z = 1e20
        max_z = 0
        for line in gro_i[1]:
            z = float(line.split()[-1])
            if z < min_z:
                min_z = z
            if z > max_z:
                max_z = z
        delta_ar.append(max_z-min_z)
    delta_ar = np.array(delta_ar); delta_ar += 0.2 # Increase the delta for each monomer by a small amount to prevent overlapping

    delta = 0
    for i, poly in enumerate(Copoly):
        if i == 0: # Head
            groFile = gro[poly][0]
            #groFile = np.append(groFile, ';Head\n')

            num = mon_len[poly][0]; delta += delta_ar[poly]

        elif i == b_Tot-1: # Tail
            edGro = editGro(gro[poly][2], num, delta, i)
            for line in edGro:
                groFile = np.append(groFile, line)
            #groFile = np.append(groFile, ';Tail\n')

        else: # Mid
            edGro = editGro(gro[poly][1], num, delta, i)
            for line in edGro:
                groFile = np.append(groFile, line)
            #groFile = np.append(groFile, ';Mid\n')

            num += mon_len[poly][1]; delta += delta_ar[poly]
    
    # Adds system name and number of atoms to the beginning, and system size to the end
    groFile = np.concatenate((['{}\n'.format('_'.join(poly_names)), '{}\n'.format(len(groFile))], groFile, ['0.0 0.0 0.0\n']))
    
    #writeFile('test.gro', groFile)
    return groFile



def editAtoms(atoms, atom_num, count):
# Edits [atoms] lines in the .top file to update the atom number, residue number, and charge (apply charge scale)
#
# Inputs: atoms = array of [atoms] lines for the polymer unit; atom_num = number of atoms to be added so the numbers are sequential across the whole file;
#         count = number of residues to be added so the numbers are sequential across the whole file

    edAtom  = []
    for i, line in enumerate(atoms):
        nr = int(line.split()[0]); res_nr = int(line.split()[2]); ch = float(line.split()[6])

        edAtom.append('{:6d}'.format(nr+atom_num) + line[6:18] + '{:6d}'.format(res_nr+count) + line[24:45] + '{: 21.15f}'.format(ch*charge_scale) + line[66:])

    return edAtom

def MakeAtoms(Copoly, atoms, mon_len):
# Assembles the [atoms] section of the .top file for the polymer
#
# Inputs: Copoly = sequence of poly unit indices, where the indices is related to poly_names; atoms = array of [atoms] chunks for [poly_names,head mid tail];
#         mon_len = array of lengths of each unit monomer from poly_names

    for i, poly in enumerate(Copoly):
        if i == 0: # Head
            atomsFile = atoms[poly][0][:2]

            edAtom = editAtoms(atoms[poly][0][2:], 0, i)
            for line in edAtom:
                atomsFile = np.append(atomsFile, line)
            #atomsFile = np.append(atomsFile, ';Head\n')

            num = mon_len[poly][0]

        elif i == b_Tot-1: # Tail
            edAtom = editAtoms(atoms[poly][2][2:], num, i)
            for line in edAtom:
                atomsFile = np.append(atomsFile, line)
            #atomsFile = np.append(atomsFile, ';Tail\n')

        else: # Mid
            edAtom = editAtoms(atoms[poly][1][2:], num, i)
            for line in edAtom:
                atomsFile = np.append(atomsFile, line)
            #atomsFile = np.append(atomsFile, ';Mid\n')

            num += mon_len[poly][1]
    
    atomsFile = np.concatenate((atomsFile, ['\n\n\n']))
    
    #writeFile('test.atoms', atomsFile)
    return atomsFile



def editBonds(bonds, atom_num):
# Edits [bonds] lines in the .top file to update the atom numbers
#
# Inputs: bonds = array of [bonds] lines for the polymer unit; atom_num = number of atoms to be added so the numbers match with the [atoms] section
    edBond  = []
    for i, line in enumerate(bonds):
        nr = np.array(line.split()[0:3], dtype = int)
        
        nr[:-1] += atom_num

        edBond.append('{:4d} {:4d} {:d}'.format(nr[0], nr[1], nr[2]) + line[11:])
    
    return edBond

def MakeBonds(Copoly, bonds, mon_len):
# Assembles the [bonds] section of the .top file for the polymer
#
# Inputs: Copoly = sequence of poly unit indices, where the indices is related to poly_names; bonds = array of [bonds] chunks for [poly_names,head mid tail];
#         mon_len = array of lengths of each unit monomer from poly_names

    for i, poly in enumerate(Copoly):
        if i == 0: # Head
            bondFile = bonds[poly][0]
            bondFile = np.append(bondFile, ';Head\n')

            # Add the Head-Mid parameters
            for line in bonds[poly][3][1:]:
                bondFile = np.append(bondFile, line)
            bondFile = np.append(bondFile, ';Head_Mid\n')

            num = mon_len[poly][0]

        elif i == b_Tot-1: # Tail
            edBond = editBonds(bonds[poly][2][1:], num)
            for line in edBond:
                bondFile = np.append(bondFile, line)
            bondFile = np.append(bondFile, ';Tail\n')

        else: # Mid
            edBond = editBonds(bonds[poly][1][1:], num)
            for line in edBond:
                bondFile = np.append(bondFile, line)
            bondFile = np.append(bondFile, ';Mid\n')
            
            if i == b_Tot-2:
                # Add the Mid-Tail parameters
                edBond = editBonds(bonds[poly][4][1:], num)
                for line in edBond:
                    bondFile = np.append(bondFile, line)
                bondFile = np.append(bondFile, ';Mid_Tail\n')
            else:
                # Add the Mid-Mid parameters
                edBond = editBonds(bonds[poly][4][1:], num)
                for line in edBond:
                    bondFile = np.append(bondFile, line)
                bondFile = np.append(bondFile, ';Mid_Mid\n')
            
            num += mon_len[poly][1]
    
    bondFile = np.concatenate((bondFile, ['\n\n\n']))
    
    #writeFile('test.bonds', bondFile)
    return bondFile



def editAngles(angles, atom_num, mon_len, poly, poly_next, check):
# Edits [angles] lines in the .top file to update the atom number
#
# Inputs: angles = array of [angles] lines for the polymer unit; atom_num = number of atoms to be added so the numbers match with the [atoms] section
#         mon_len = array of lengths of each unit monomer from poly_names; poly = identity of current monomer, poly_next = identity of next monomer, check = 0 if Mid-Mid, 1 if Mid-Tail

    edAngle  = []
    for i, line in enumerate(angles):
        nr = np.array(line.split()[0:4], dtype = int)

        if len(nr[nr > 2*mon_len[poly][1]]) > 0: # Deals with Mid-Tail
            if check == 0: # Skip Mid-Tail parameters (only those that interact with the tail cap (e.g. Hydrogen capping an ethylene backhone)) if Mid-Mid Section
                continue
            elif check == 1: # Edit Mid-Tail parameters correctly if Mid-Tail section
                nr[nr > 2*mon_len[poly][1]] -= (mon_len[poly][1] - mon_len[poly_next][1])

        nr[:-1] += atom_num

        edAngle.append('{:4d} {:4d} {:4d} {:d}'.format(nr[0], nr[1], nr[2], nr[3]) + line[16:])
    
    return edAngle

def editAnglesMMM(angles, atom_num, mon_len, poly_prev, poly, index): #NOTE: This function is currently untested, however it functions exactly like the following MMM functions, so it should work fine
# Edits Head-Mid-Mid/Mid-Mid-Mid/Mid-Mid-Tail [angles] lines in the .top file to update the atom number
#
# Inputs = angles = array of [angles] lines for the polymer unit; atom_num = number of atoms to be added so the numbers match with the [atoms] section
#          mon_len = array of lengths of each unit monomer from poly_names; poly_prev = identity of previous monomer, poly = identity of current monomer, index = index of Copoly
    
    edAngle  = []
    for i, line in enumerate(angles):

        nr = np.array(line.split()[0:4], dtype = int); nbr = nr[-1]
        nr[:-1] -= mon_len[poly][1]

        if index == 1: # Deals with Head-Mid-Mid
            nr[nr < 1] += (mon_len[poly][1] - mon_len[poly_prev][0])
        else: # Deals with Mid-Mid-Mid
            nr[nr < 1] += (mon_len[poly][1] - mon_len[poly_prev][1])

        nr[nr > 0] += (mon_len[poly][1] - mon_len[poly][0])
        nr[:-1] += atom_num

        edAngle.append('{:4d} {:4d} {:4d} {:d}'.format(nr[0], nr[1], nr[2], nbr) + line[16:])
    
    return edAngle

def MakeAngles(Copoly, angles, mon_len):
# Assembles the [angles] section of the .top file for the polymer
#
# Inputs: Copoly = sequence of poly unit indices, where the indices is related to poly_names; angles = array of [angles] chunks for [poly_names,head mid tail];
#         mon_len = array of lengths of each unit monomer from poly_names

    for i, poly in enumerate(Copoly):
        if i == 0: # Head
            anglesFile = angles[poly][0]
            anglesFile = np.append(anglesFile, ';Head\n')

            # Add the Head-Mid parameters
            for line in angles[poly][3][1:]:
                anglesFile = np.append(anglesFile, line)
            anglesFile = np.append(anglesFile, ';Head_Mid\n')
            
            num = mon_len[poly][0]
        elif i == b_Tot-1: # Tail
            edAngle = editAngles(angles[poly][2][1:], num, mon_len, poly, poly, 0)
            for line in edAngle:
                anglesFile = np.append(anglesFile, line)
            anglesFile = np.append(anglesFile, ';Tail\n')
        else: # Mid
            edAngle = editAngles(angles[poly][1][1:], num, mon_len, poly, Copoly[i+1], 0)
            for line in edAngle:
                anglesFile = np.append(anglesFile, line)
            anglesFile = np.append(anglesFile, ';Mid\n')

            if i == b_Tot-2:
                # Add the Mid-Tail parameters
                edAngle = editAngles(angles[poly][4][1:], num, mon_len, poly, Copoly[i+1], 1)
                for line in edAngle:
                    anglesFile = np.append(anglesFile, line)
                anglesFile = np.append(anglesFile, ';Mid_Tail\n')
            else:
                # Add the Mid-Mid parameters
                edAngle = editAngles(angles[poly][4][1:], num, mon_len, poly, Copoly[i+1], 0)
                for line in edAngle:
                    anglesFile = np.append(anglesFile, line)
                anglesFile = np.append(anglesFile, ';Mid_Mid\n')
            
            # Add the Head-Mid-Mid/Mid-Mid-Mid/Mid-Mid-Tail parameters
            edAngle = editAnglesMMM(angles[poly][5][1:], num, mon_len, Copoly[i-1], poly, i)
            for line in edAngle:
                anglesFile = np.append(anglesFile, line)
            anglesFile = np.append(anglesFile, ';Mid_Mid_Mid\n')
            
            num += mon_len[poly][1]
    
    anglesFile = np.concatenate((anglesFile, ['\n\n\n']))
    
    #writeFile('test.angles', anglesFile)
    return anglesFile



def editDihedrals(dihedrals, atom_num, mon_len, poly, poly_next, check):
# Edits [dihedrals] lines in the .top file to update the atom number
#
# Inputs: dihedrals = array of [dihedrals] lines for the polymer unit; atom_num = number of atoms to be added so the numbers match with the [atoms] section
#         mon_len = array of lengths of each unit monomer from poly_names; poly = identity of current monomer, poly_next = identity of next monomer, check = 0 if Mid-Mid, 1 if Mid-Tail

    edDihedral  = []
    for i, line in enumerate(dihedrals):
        nr = np.array(line.split()[0:5], dtype = int)

        if len(nr[nr > 2*mon_len[poly][1]]) > 0: # Deals with Mid-Tail
            if check == 0: # Skip Mid-Tail parameters (only those that interact with the tail cap (e.g. Hydrogen capping an ethylene backhone)) if Mid-Mid Section
                continue
            elif check == 1: # Edit Mid-Tail parameters correctly if Mid-Tail section
                nr[nr > 2*mon_len[poly][1]] -= (mon_len[poly][1] - mon_len[poly_next][1])

        nr[:-1] += atom_num

        edDihedral.append('{:4d} {:4d} {:4d} {:4d} {:d}'.format(nr[0], nr[1], nr[2], nr[3], nr[4]) + line[21:])
    
    return edDihedral

def editDihedralsMMM(dihedrals, atom_num, mon_len, poly_prev, poly, index):
# Edits Head-Mid-Mid/Mid-Mid-Mid/Mid-Mid-Tail [dihedrals] lines in the .top file to update the atom number
#
# Inputs = dihedrals = array of [dihedrals] lines for the polymer unit; atom_num = number of atoms to be added so the numbers match with the [atoms] section
#          mon_len = array of lengths of each unit monomer from poly_names; poly_prev = identity of previous monomer, poly = identity of current monomer, index = index of Copoly

    edDihedral  = []
    for i, line in enumerate(dihedrals):

        nr = np.array(line.split()[0:5], dtype = int); nbr = nr[-1]
        nr[:-1] -= mon_len[poly][1]

        if index == 1: # Deals with Head-Mid-Mid
            nr[nr < 1] += (mon_len[poly][1] - mon_len[poly_prev][0])
        else: # Deals with Mid-Mid-Mid
            nr[nr < 1] += (mon_len[poly][1] - mon_len[poly_prev][1])

        nr[nr > 0] += (mon_len[poly][1] - mon_len[poly][0])
        nr[:-1] += atom_num

        edDihedral.append('{:4d} {:4d} {:4d} {:4d} {:d}'.format(nr[0], nr[1], nr[2], nr[3], nbr) + line[21:])
    
    return edDihedral

def MakeDihedrals(Copoly, dihedrals, mon_len):
# Assembles the [dihedrals] section of the .top file for the polymer
#
# Inputs: Copoly = sequence of poly unit indices, where the indices is related to poly_names; dihedrals = array of [atoms] chunks for [poly_names,head mid tail];
#         mon_len = array of lengths of each unit monomer from poly_names

    for i, poly in enumerate(Copoly):
        if i == 0: # Head
            dihedralsFile = dihedrals[poly][0]
            dihedralsFile = np.append(dihedralsFile, ';Head\n')

            # Add the Head-Mid parameters
            for line in dihedrals[poly][3][1:]:
                dihedralsFile = np.append(dihedralsFile, line)
            dihedralsFile = np.append(dihedralsFile, ';Head_Mid\n')
            
            num = mon_len[poly][0]

        elif i == b_Tot-1: # Tail
            edDihedrals = editDihedrals(dihedrals[poly][2][1:], num, mon_len, poly, poly, 0)
            for line in edDihedrals:
                dihedralsFile = np.append(dihedralsFile, line)
            dihedralsFile = np.append(dihedralsFile, ';Tail\n')

        else: # Mid
            edDihedrals = editDihedrals(dihedrals[poly][1][1:], num, mon_len, poly, Copoly[i+1], 0)
            for line in edDihedrals:
                dihedralsFile = np.append(dihedralsFile, line)
            dihedralsFile = np.append(dihedralsFile, ';Mid\n')
            
            if i == b_Tot-2:
                # Add the Mid-Tail parameters
                edDihedrals = editDihedrals(dihedrals[poly][4][1:], num, mon_len, poly, Copoly[i+1], 1)
                for line in edDihedrals:
                    dihedralsFile = np.append(dihedralsFile, line)
                dihedralsFile = np.append(dihedralsFile, ';Mid_Tail\n')
            else:
                # Add the Mid-Mid parameters
                edDihedrals = editDihedrals(dihedrals[poly][4][1:], num, mon_len, poly, Copoly[i+1], 0)
                for line in edDihedrals:
                    dihedralsFile = np.append(dihedralsFile, line)
                dihedralsFile = np.append(dihedralsFile, ';Mid_Mid\n')
            
            # Add the Head-Mid-Mid/Mid-Mid-Mid/Mid-Mid-Tail parameters
            edDihedrals = editDihedralsMMM(dihedrals[poly][5][1:], num, mon_len, Copoly[i-1], poly, i)
            for line in edDihedrals:
                dihedralsFile = np.append(dihedralsFile, line)
            dihedralsFile = np.append(dihedralsFile, ';Mid_Mid_Mid\n')
            
            num += mon_len[poly][1]
    
    dihedralsFile = np.concatenate((dihedralsFile, ['\n\n\n']))
    
    #writeFile('test.dihedrals', dihedralsFile)
    return dihedralsFile



def editImpropers(impropers, atom_num, mon_len, poly, poly_next, check):
# Edits [improper] lines in the .top file to update the atom number
#
# Inputs: impropers = array of [improper] lines for the polymer unit; atom_num = number of atoms to be added so the numbers match with the [atoms] section
#         mon_len = array of lengths of each unit monomer from poly_names; poly = identity of current monomer, poly_next = identity of next monomer, check = 0 if Mid-Mid, 1 if Mid-Tail

    edImproper  = []
    for i, line in enumerate(impropers):
        nr = np.array(line.split()[0:5], dtype = int)
    
        if len(nr[nr > 2*mon_len[poly][1]]) > 0: # Deals with Mid-Tail
            if check == 0: # Skip Mid-Tail parameters (only those that interact with the tail cap (e.g. Hydrogen capping an ethylene backhone)) if Mid-Mid Section
                continue
            elif check == 1: # Edit Mid-Tail parameters correctly if Mid-Tail section
                nr[nr > 2*mon_len[poly][1]] -= (mon_len[poly][1] - mon_len[poly_next][1])

        nr[:-1] += atom_num

        edImproper.append('{:4d} {:4d} {:4d} {:4d} {:d}'.format(nr[0], nr[1], nr[2], nr[3], nr[4]) + line[21:])
    
    return edImproper

def editImpropersMMM(impropers, atom_num, mon_len, poly_prev, poly, index):
# Edits Head-Mid-Mid/Mid-Mid-Mid/Mid-Mid-Tail [impropers] lines in the .top file to update the atom number
#
# Inputs = impropers = array of [impropers] lines for the polymer unit; atom_num = number of atoms to be added so the numbers match with the [atoms] section
#          mon_len = array of lengths of each unit monomer from poly_names; poly_prev = identity of previous monomer, poly = identity of current monomer, index = index of Copoly

    edImproper  = []
    for i, line in enumerate(impropers):

        nr = np.array(line.split()[0:5], dtype = int); nbr = nr[-1]
        nr[:-1] -= mon_len[poly][1]

        if index == 1: # Deals with Head-Mid-Mid
            nr[nr < 1] += (mon_len[poly][1] - mon_len[poly_prev][0])
        else: # Deals with Mid-Mid-Mid
            nr[nr < 1] += (mon_len[poly][1] - mon_len[poly_prev][1])

        nr[nr > 0] += (mon_len[poly][1] - mon_len[poly][0])
        nr[:-1] += atom_num

        edImproper.append('{:4d} {:4d} {:4d} {:4d} {:d}'.format(nr[0], nr[1], nr[2], nr[3], nbr) + line[21:])
    
    return edImproper

def MakeImpropers(Copoly, impropers, mon_len):
# Assembles the [impropers] section of the .top file for the polymer
#
# Inputs: Copoly = sequence of poly unit indices, where the indices is related to poly_names; impropers = array of [atoms] chunks for [poly_names,head mid tail];
#         mon_len = array of lengths of each unit monomer from poly_names

    for i, poly in enumerate(Copoly):
        if i == 0: # Head
            impropersFile = impropers[poly][0]
            impropersFile = np.append(impropersFile, ';Head\n')

            # Add the Head-Mid parameters
            for line in impropers[poly][3][1:]:
                impropersFile = np.append(impropersFile, line)
            impropersFile = np.append(impropersFile, ';Head_Mid\n')
            
            num = mon_len[poly][0]

        elif i == b_Tot-1: # Tail
            edImpropers = editImpropers(impropers[poly][2][1:], num, mon_len, poly, poly, 0)
            for line in edImpropers:
                impropersFile = np.append(impropersFile, line)
            impropersFile = np.append(impropersFile, ';Tail\n')

        else: # Mid
            edImpropers = editImpropers(impropers[poly][1][1:], num, mon_len, poly, Copoly[i+1], 0)
            for line in edImpropers:
                impropersFile = np.append(impropersFile, line)
            impropersFile = np.append(impropersFile, ';Mid\n')

            if i == b_Tot-2:
                # Add the Mid-Tail parameters
                edImpropers = editImpropers(impropers[poly][4][1:], num, mon_len, poly, Copoly[i+1], 1)
                for line in edImpropers:
                    impropersFile = np.append(impropersFile, line)
                impropersFile = np.append(impropersFile, ';Mid_Tail\n')
            else:
                # Add the Mid-Mid parameters
                edImpropers = editImpropers(impropers[poly][4][1:], num, mon_len, poly, Copoly[i+1], 0)
                for line in edImpropers:
                    impropersFile = np.append(impropersFile, line)
                impropersFile = np.append(impropersFile, ';Mid_Mid\n')
            
            # Add the Head-Mid-Mid/Mid-Mid-Mid/Mid-Mid-Tail parameters
            edImpropers = editImpropersMMM(impropers[poly][5][1:], num, mon_len, Copoly[i-1], poly, i)
            for line in edImpropers:
                impropersFile = np.append(impropersFile, line)
            impropersFile = np.append(impropersFile, ';Mid_Mid_Mid\n')
            
            num += mon_len[poly][1]
    
    impropersFile = np.concatenate((impropersFile, ['\n\n\n']))
    
    #writeFile('test.impropers', impropersFile)
    return impropersFile



def editPairs(pairs, atom_num, mon_len, poly, poly_next, check):
# Edits [pairs] lines in the .top file to update the atom number
#
# Inputs: pairs = array of [pairs] lines for the polymer unit; atom_num = number of atoms to be added so the numbers match with the [atoms] section
#         mon_len = array of lengths of each unit monomer from poly_names; poly = identity of current monomer, poly_next = identity of next monomer, check = 0 if Mid-Mid, 1 if Mid-Tail

    edPair  = []
    for i, line in enumerate(pairs):
        nr = np.array(line.split()[0:3], dtype = int)

        if len(nr[nr > 2*mon_len[poly][1]]) > 0: # Deals with Mid-Tail
            if check == 0: # Skip Mid-Tail parameters (only those that interact with the tail cap (e.g. Hydrogen capping an ethylene backhone)) if Mid-Mid Section
                continue
            elif check == 1: # Edit Mid-Tail parameters correctly if Mid-Tail section
                nr[nr > 2*mon_len[poly][1]] -= (mon_len[poly][1] - mon_len[poly_next][1])

        nr[:-1] += atom_num


        edPair.append('{:4d} {:4d} {:d}'.format(nr[0], nr[1], nr[2]) + line[11:])
    
    return edPair

def editPairsMMM(pairs, atom_num, mon_len, poly_prev, poly, index):
# Edits Head-Mid-Mid/Mid-Mid-Mid/Mid-Mid-Tail [pairs] lines in the .top file to update the atom number
#
# Inputs = pairs = array of [pairs] lines for the polymer unit; atom_num = number of atoms to be added so the numbers match with the [atoms] section
#          mon_len = array of lengths of each unit monomer from poly_names; poly_prev = identity of previous monomer, poly = identity of current monomer, index = index of Copoly

    edPair  = []
    for i, line in enumerate(pairs):

        nr = np.array(line.split()[0:3], dtype = int); nbr = nr[-1]
        nr[:-1] -= mon_len[poly][1]

        if index == 1: # Deals with Head-Mid-Mid
            nr[nr < 1] += (mon_len[poly][1] - mon_len[poly_prev][0])
        else: # Deals with Mid-Mid-Mid
            nr[nr < 1] += (mon_len[poly][1] - mon_len[poly_prev][1])

        nr[nr > 0] += (mon_len[poly][1] - mon_len[poly][0])
        nr[:-1] += atom_num

        edPair.append('{:4d} {:4d} {:d}'.format(nr[0], nr[1], nbr) + line[11:])
    
    return edPair

def MakePairs(Copoly, pairs, mon_len):
# Assembles the [pairs] section of the .top file for the polymer
#
# Inputs: Copoly = sequence of poly unit indices, where the indices is related to poly_names; pairs = array of [atoms] chunks for [poly_names,head mid tail];
#         mon_len = array of lengths of each unit monomer from poly_names

    for i, poly in enumerate(Copoly):
        if i == 0:
            pairsFile = pairs[poly][0]
            pairsFile = np.append(pairsFile, ';Head\n')

            # Add the Head-Mid parameters
            for line in pairs[poly][3][1:]:
                pairsFile = np.append(pairsFile, line)
            pairsFile = np.append(pairsFile, ';Head_Mid\n')
            
            num = mon_len[poly][0]

        elif i == b_Tot-1:
            edPairs = editPairs(pairs[poly][2][1:], num, mon_len, poly, poly, 0)
            for line in edPairs:
                pairsFile = np.append(pairsFile, line)
            pairsFile = np.append(pairsFile, ';Tail\n')

        else:
            edPairs = editPairs(pairs[poly][1][1:], num, mon_len, poly, Copoly[i+1], 0)
            for line in edPairs:
                pairsFile = np.append(pairsFile, line)
            pairsFile = np.append(pairsFile, ';Mid\n')
            
            if i == b_Tot-2:
                # Add the Mid-Tail parameters
                edPairs = editPairs(pairs[poly][4][1:], num, mon_len, poly, Copoly[i+1], 1)
                for line in edPairs:
                    pairsFile = np.append(pairsFile, line)
                pairsFile = np.append(pairsFile, ';Mid_Tail\n')
            else:
                # Add the Mid-Mid parameters
                edPairs = editPairs(pairs[poly][4][1:], num, mon_len, poly, Copoly[i+1], 0)
                for line in edPairs:
                    pairsFile = np.append(pairsFile, line)
                pairsFile = np.append(pairsFile, ';Mid_Mid\n')
            
            # Add the Head-Mid-Mid/Mid-Mid-Mid/Mid-Mid-Tail parameters
            edPairs = editPairsMMM(pairs[poly][5][1:], num, mon_len, Copoly[i-1], poly, i)
            for line in edPairs:
                pairsFile = np.append(pairsFile, line)
            pairsFile = np.append(pairsFile, ';Mid_Mid_Mid\n')
            
            num += mon_len[poly][1]
    
    pairsFile = np.concatenate((pairsFile, ['\n\n\n']))
    
    #writeFile('test.pairs', pairsFile)
    return pairsFile



def main():

    if poly_type == 0:
        Copoly = PolySeq_Rand()
    elif poly_type == 1:
        Copoly = PolySeq_Block()

    mon_len, mon_nums, head, gro, atoms, bonds, angles, dihedrals, impropers, pairs, tail = PolyRead()

    groFile = MakeGro(Copoly, gro, mon_len)

    atomsFile = MakeAtoms(Copoly, atoms, mon_len)

    bondsFile = MakeBonds(Copoly, bonds, mon_len)

    anglesFile = MakeAngles(Copoly, angles, mon_len)

    dihedralsFile = MakeDihedrals(Copoly, dihedrals, mon_len)

    impropersFile = MakeImpropers(Copoly, impropers, mon_len)

    pairsFile = MakePairs(Copoly, pairs, mon_len)

    # Add the text in includes.txt after the [pairs] section and before the tail
    # Generally used to add water and ion itp file paths
    if os.path.isfile('header.txt') and os.path.isfile('includes.txt'):
        with open('header.txt', 'r') as f:
            header = f.readlines()
        with open('includes.txt', 'r') as f:
            includes = f.readlines()

        topFile = np.concatenate((header, atomsFile, bondsFile, anglesFile, dihedralsFile, impropersFile, pairsFile, includes, tail[0], ['\n']))

    else:
        topFile = np.concatenate((head[0], atomsFile, bondsFile, anglesFile, dihedralsFile, impropersFile, pairsFile, tail[0], ['\n']))


    if not os.path.isdir('_'.join(poly_names)+'/'+folder+'/sample_'+sample):
        os.makedirs('_'.join(poly_names)+'/'+folder+'/sample_'+sample)

    writeFile('_'.join(poly_names)+'/'+folder+'/sample_'+sample+'/polymer.gro', groFile)
    writeFile('_'.join(poly_names)+'/'+folder+'/sample_'+sample+'/topol.top', topFile)

if __name__ == "__main__":
    main()
