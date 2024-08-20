import numpy as np
import h5py
import os

from sys import argv
script, poly_name = argv
# folder = folder containing trimer.top file

os.chdir(poly_name)
names = ['h', 'm', 't', 'hm', 'mt', 'hmt']


# Read in .gro files
# NOTE: .gro files should be in the format made from gmx pdb2gmx
with open(poly_name+'h.gro') as f:
    gro_h = f.readlines()

    h_len = int(gro_h[1]) # Number of Atoms
    gro_h = np.array(gro_h); gro_h = gro_h[2:-1]

with open(poly_name+'m.gro') as f:
    gro_m = f.readlines()

    m_len = int(gro_m[1]) # Number of Atoms
    gro_m = np.array(gro_m); gro_m = gro_m[2:-1]

with open(poly_name+'t.gro') as f:
    gro_t = f.readlines()

    t_len = int(gro_t[1]) # Number of Atoms
    gro_t = np.array(gro_t); gro_t = gro_t[2:-1]
gro = [gro_h, gro_m, gro_t]
#print(gro_h)
#print(gro_m)
#print(gro_t)



# Read in trimer.top and divide into sections based on headers below
# NOTE: .top file should be in the format made by my modified MKTOP.pl script
with open(poly_name+'_tri.top') as f:
    top = f.readlines()

    top = np.array(top); top = top[11:]

    head = top[:np.where(top == '[ atoms ]\n')[0][0]]
    atoms = top[np.where(top == '[ atoms ]\n')[0][0]:np.where(top == '[ bonds ]\n')[0][0]-4]
    bonds = top[np.where(top == '[ bonds ]\n')[0][0]:np.where(top == '[ angles ]\n')[0][0]-2]
    angles = top[np.where(top == '[ angles ]\n')[0][0]:np.where(top == '[ dihedrals ]\n')[0][0]-2]
    dihedrals = top[np.where(top == '[ dihedrals ]\n')[0][0]:np.where(top == '[ dihedrals ]\n')[0][1]-2]
    improper = top[np.where(top == '[ dihedrals ]\n')[0][1]:np.where(top == '[ pairs ]\n')[0][0]-2]
    pairs = top[np.where(top == '[ pairs ]\n')[0][0]:np.where(top == '[ system ]\n')[0][0]-2]
    tail = top[np.where(top == '[ system ]\n')[0][0]:]

# Retreieve the number of atoms in the monomer and the coordinate of the last atom in each unit head, mid, tail
mon_len = [h_len, m_len, t_len]                    # Number of Atoms [Head, Mid, Tail]
mon_nums = [h_len, h_len+m_len, h_len+m_len+t_len] # Last Atom Number in trimer.top [Head, Mid, Tail]
print(mon_len, mon_nums)



# Create an atoms section for the head, mid, tail units where numbering is reset as if it is the .top file for the single unit
# NOTE: cgnr is the Charge Group Number. This is no longer used for OPLSAA in GROMACS.
max_cgnr = 1
atoms_ar = [atoms[:2], atoms[:2], atoms[:2]] # [Head, Mid, Tail]
for line in atoms[2:]:
    nr = int(line.split()[0]); cgnr = int(line.split()[5])
    if nr <= mon_nums[0]:
        line = line[:24] + '{:>8s}'.format(line.split()[3]) + line[32:]
        atoms_ar[0] = np.append(atoms_ar[0], line)

        if cgnr > max_cgnr:
            max_cgnr = cgnr
    elif nr <= mon_nums[1]:
        line = '{:6d}'.format(nr-mon_nums[0]) + line[6:38] + '{:7d}'.format(cgnr-max_cgnr) + line[45:]
        atoms_ar[1] = np.append(atoms_ar[1], line)
    else:
        line = line[:24] + '{:>8s}'.format(line.split()[3]) + line[32:]
        line = '{:6d}'.format(nr-mon_nums[1]) + line[6:38] + '{:7d}'.format(cgnr-2*max_cgnr) + line[45:]
        atoms_ar[2] = np.append(atoms_ar[2], line)



# Create a bonds section for the head, mid, tail units where numbering is reset as if it is the .top file for the single unit
# Bonds between units are found and saved separately from each unit
bonds_ar = [bonds[:1], bonds[:1], bonds[:1], bonds[:1], bonds[:1], bonds[:1]] # [Head, Mid, Tail, Head-Mid, Mid-Tail, Head-Mid-Tail]
for line in bonds[1:]:
    nrs = np.array(line.split()[:2], dtype = int)
    if np.all(nrs <= mon_nums[0]):
        bonds_ar[0] = np.append(bonds_ar[0], line)
    elif np.all(nrs > mon_nums[0]) and np.all(nrs <= mon_nums[1]):
        line = '{:4d} {:4d}'.format(nrs[0]-mon_nums[0], nrs[1]-mon_nums[0]) + line[9:]
        bonds_ar[1] = np.append(bonds_ar[1], line)
    elif np.all(nrs > mon_nums[1]) and np.all(nrs <= mon_nums[2]):
        line = '{:4d} {:4d}'.format(nrs[0]-mon_nums[1], nrs[1]-mon_nums[1]) + line[9:]
        bonds_ar[2] = np.append(bonds_ar[2], line)
    elif np.all(nrs <= mon_nums[1]):
        bonds_ar[3] = np.append(bonds_ar[3], line)
    elif not np.all(nrs > mon_nums[0]):
        bonds_ar[5] = np.append(bonds_ar[5], line)
    else:
        line = '{:4d} {:4d}'.format(nrs[0]-mon_nums[0], nrs[1]-mon_nums[0]) + line[9:]
        bonds_ar[4] = np.append(bonds_ar[4], line)
#print(bonds_ar[3])
#print(bonds_ar[4])
#print(bonds_ar[5])



# Create an angles section for the head, mid, tail units where numbering is reset as if it is the .top file for the single unit
# Angles between units are found and saved separately from each unit
angles_ar = [angles[:1], angles[:1], angles[:1], angles[:1], angles[:1], angles[:1]] # [Head, Mid, Tail, Head-Mid, Mid-Tail, Head-Mid-Tail]
for line in angles[1:]:
    nrs = np.array(line.split()[:3], dtype = int)
    if np.all(nrs <= mon_nums[0]):
        angles_ar[0] = np.append(angles_ar[0], line)
    elif np.all(nrs > mon_nums[0]) and np.all(nrs <= mon_nums[1]):
        line = '{:4d} {:4d} {:4d}'.format(nrs[0]-mon_nums[0], nrs[1]-mon_nums[0], nrs[2]-mon_nums[0]) + line[14:]
        angles_ar[1] = np.append(angles_ar[1], line)
    elif np.all(nrs > mon_nums[1]) and np.all(nrs <= mon_nums[2]):
        line = '{:4d} {:4d} {:4d}'.format(nrs[0]-mon_nums[1], nrs[1]-mon_nums[1], nrs[2]-mon_nums[1]) + line[14:]
        angles_ar[2] = np.append(angles_ar[2], line)
    elif np.all(nrs <= mon_nums[1]):
        angles_ar[3] = np.append(angles_ar[3], line)
    elif not np.all(nrs > mon_nums[0]):
        angles_ar[5] = np.append(angles_ar[5], line)
    else:
        line = '{:4d} {:4d} {:4d}'.format(nrs[0]-mon_nums[0], nrs[1]-mon_nums[0], nrs[2]-mon_nums[0]) + line[14:]
        angles_ar[4] = np.append(angles_ar[4], line)
#print(angles_ar[3])
#print(angles_ar[4])
#print(angles_ar[5])



# Create a dihedrals section for the head, mid, tail units where numbering is reset as if it is the .top file for the single unit
# Dihedrals between units are found and saved separately from each unit
dihedrals_ar = [dihedrals[:1], dihedrals[:1], dihedrals[:1], dihedrals[:1], dihedrals[:1], dihedrals[:1]] # [Head, Mid, Tail, Head-Mid, Mid-Tail, Head-Mid-Tail]
for line in dihedrals[1:]:
    nrs = np.array(line.split()[:4], dtype = int)
    if np.all(nrs <= mon_nums[0]):
        dihedrals_ar[0] = np.append(dihedrals_ar[0], line)
    elif np.all(nrs > mon_nums[0]) and np.all(nrs <= mon_nums[1]):
        line = '{:4d} {:4d} {:4d} {:4d}'.format(nrs[0]-mon_nums[0], nrs[1]-mon_nums[0], nrs[2]-mon_nums[0], nrs[3]-mon_nums[0]) + line[19:]
        dihedrals_ar[1] = np.append(dihedrals_ar[1], line)
    elif np.all(nrs > mon_nums[1]) and np.all(nrs <= mon_nums[2]):
        line = '{:4d} {:4d} {:4d} {:4d}'.format(nrs[0]-mon_nums[1], nrs[1]-mon_nums[1], nrs[2]-mon_nums[1], nrs[3]-mon_nums[1]) + line[19:]
        dihedrals_ar[2] = np.append(dihedrals_ar[2], line)
    elif np.all(nrs <= mon_nums[1]):
        dihedrals_ar[3] = np.append(dihedrals_ar[3], line)
    elif not np.all(nrs > mon_nums[0]):
        dihedrals_ar[5] = np.append(dihedrals_ar[5], line)
    else:
        line = '{:4d} {:4d} {:4d} {:4d}'.format(nrs[0]-mon_nums[0], nrs[1]-mon_nums[0], nrs[2]-mon_nums[0], nrs[3]-mon_nums[0]) + line[19:]
        dihedrals_ar[4] = np.append(dihedrals_ar[4], line)
#print(dihedrals_ar[3])
#print(dihedrals_ar[4])
#print(dihedrals_ar[5])



# Create an improper dihedrals section for the head, mid, tail units where numbering is reset as if it is the .top file for the single unit
# Improper Dihedrals between units are found and saved separately from each unit
improper_ar = [improper[:1], improper[:1], improper[:1], improper[:1], improper[:1], improper[:1]] # [Head, Mid, Tail, Head-Mid, Mid-Tail, Head-Mid-Tail]
for line in improper[1:]:
    nrs = np.array(line.split()[:4], dtype = int)
    if np.all(nrs <= mon_nums[0]):
        improper_ar[0] = np.append(improper_ar[0], line)
    elif np.all(nrs > mon_nums[0]) and np.all(nrs <= mon_nums[1]):
        line = '{:4d} {:4d} {:4d} {:4d}'.format(nrs[0]-mon_nums[0], nrs[1]-mon_nums[0], nrs[2]-mon_nums[0], nrs[3]-mon_nums[0]) + line[19:]
        improper_ar[1] = np.append(improper_ar[1], line)
    elif np.all(nrs > mon_nums[1]) and np.all(nrs <= mon_nums[2]):
        line = '{:4d} {:4d} {:4d} {:4d}'.format(nrs[0]-mon_nums[1], nrs[1]-mon_nums[1], nrs[2]-mon_nums[1], nrs[3]-mon_nums[1]) + line[19:]
        improper_ar[2] = np.append(improper_ar[2], line)
    elif np.all(nrs <= mon_nums[1]):
        improper_ar[3] = np.append(improper_ar[3], line)
    elif not np.all(nrs > mon_nums[0]):
        improper_ar[5] = np.append(improper_ar[5], line)
    else:
        line = '{:4d} {:4d} {:4d} {:4d}'.format(nrs[0]-mon_nums[0], nrs[1]-mon_nums[0], nrs[2]-mon_nums[0], nrs[3]-mon_nums[0]) + line[19:]
        improper_ar[4] = np.append(improper_ar[4], line)
#print(improper_ar[3])
#print(improper_ar[4])
#print(improper_ar[5])



# Create a pairs section for the head, mid, tail units where numbering is reset as if it is the .top file for the single unit
# Pairs between units are found and saved separately from each unit
pairs_ar = [pairs[:1], pairs[:1], pairs[:1], pairs[:1], pairs[:1], pairs[:1]] # [Head, Mid, Tail, Head-Mid, Mid-Tail, Head-Mid-Tail]
for line in pairs[1:]:
    nrs = np.array(line.split()[:2], dtype = int)
    if np.all(nrs <= mon_nums[0]):
        pairs_ar[0] = np.append(pairs_ar[0], line)
    elif np.all(nrs > mon_nums[0]) and np.all(nrs <= mon_nums[1]):
        line = '{:4d} {:4d}'.format(nrs[0]-mon_nums[0], nrs[1]-mon_nums[0]) + line[9:]
        pairs_ar[1] = np.append(pairs_ar[1], line)
    elif np.all(nrs > mon_nums[1]) and np.all(nrs <= mon_nums[2]):
        line = '{:4d} {:4d}'.format(nrs[0]-mon_nums[1], nrs[1]-mon_nums[1]) + line[9:]
        pairs_ar[2] = np.append(pairs_ar[2], line)
    elif np.all(nrs <= mon_nums[1]):
        pairs_ar[3] = np.append(pairs_ar[3], line)
    elif not np.all(nrs > mon_nums[0]):
        pairs_ar[5] = np.append(pairs_ar[5], line)
    else:
        line = '{:4d} {:4d}'.format(nrs[0]-mon_nums[0], nrs[1]-mon_nums[0]) + line[9:]
        pairs_ar[4] = np.append(pairs_ar[4], line)
#print(pairs_ar[3])
#print(pairs_ar[4])
#print(pairs_ar[5])



# Save all .gro and .top data to a h5py file for use in the poly_build.py code
with h5py.File('top.hdf5','w') as f:
    dset1 = f.create_dataset("mon_len", data=mon_len)
    dset2 = f.create_dataset("mon_nums", data=mon_nums)
    dset3 = f.create_dataset("head", data=np.array(head, dtype=object))
    for i, n in enumerate(names):
        if i < 3:
            dset4 = f.create_dataset(n+"_atoms", data=np.array(atoms_ar[i], dtype=object))
        dset5 = f.create_dataset(n+"_bonds", data=np.array(bonds_ar[i], dtype=object))
        dset6 = f.create_dataset(n+"_angles", data=np.array(angles_ar[i], dtype=object))
        dset7 = f.create_dataset(n+"_dihedrals", data=np.array(dihedrals_ar[i], dtype=object))
        dset8 = f.create_dataset(n+"_improper", data=np.array(improper_ar[i], dtype=object))
        dset9 = f.create_dataset(n+"_pairs", data=np.array(pairs_ar[i], dtype=object))

        if i < 3:
            dset10 = f.create_dataset(n+"_gro", data=np.array(gro[i], dtype=object))
    dset11 = f.create_dataset("tail", data=np.array(tail, dtype=object))
