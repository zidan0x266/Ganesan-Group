# This script generates a single topology file containing all of the monomers. This topology file will be used for crosslinking the monomers

import builder_io
import copy
import numpy as np
import os

# If new TOP exists, delete it
if os.path.isfile('polyamide.top'):
    os.remove('polyamide.top')

# Read in molecule TOPS + starter TOP
TOP = builder_io.GroTopFile('MPD.top'); TOP.read()
MPD = builder_io.GroTopFile('MPD.top'); MPD.read()
TMC = builder_io.GroTopFile('TMC.top'); TMC.read()

# Set basic TOP parameters
TOP.system_name = 'Polyamide'
TOP.moleculetype = [{'name': 'MOL', 'nrexcl': '3'}]
TOP.molecules = [{'name': 'MOL', 'mol': '1'}]

# Number of molecules; number of atoms per molecule, name of crosslinker atoms
nMPD = 690; lMPD = len(MPD.atoms)
nTMC = 460; lTMC = len(TMC.atoms)



### Atoms Section ###
# Molecule 1
n = len(TOP.atoms) + 1
for i in range(2, nMPD + 1):
    for j in range(1, len(MPD.atoms) + 1):
        # Add new atom, update atom and chain index
        TOP.atoms[n] = copy.copy(MPD.atoms[j])
        TOP.atoms[n].atom_id = n
        TOP.atoms[n].chain_idx = i
        n += 1
#Molecule 2
for i in range(1, nTMC + 1):
    for j in range(1, len(TMC.atoms) + 1):
        # Add new atom, update atom and chain index
        TOP.atoms[n] = copy.copy(TMC.atoms[j])
        TOP.atoms[n].atom_id = n
        TOP.atoms[n].chain_idx = nMPD + i
        n += 1



### Bonds Section ###
# Molecule 1
for i in range(2, nMPD + 1):
    bMPD = copy.copy(MPD.bonds)
    for j in bMPD:
        # Add new bond with updated atom indexes
        TOP.bonds[tuple(np.array(j) + lMPD*(i-1))] = copy.copy(bMPD[j])
# Molecule 2
for i in range(1, nTMC + 1):
    bTMC = copy.copy(TMC.bonds)
    for j in bTMC:
        # Add new bond with updated atom indexes
        TOP.bonds[tuple(np.array(j) + lTMC*(i-1) + lMPD*nMPD)] = copy.copy(bTMC[j])



### Angles Section ###
# Molecule 1
for i in range(2, nMPD + 1):
    aMPD = copy.copy(MPD.angles)
    for j in aMPD:
        # Add new angle with updated atom indexes
        TOP.angles[tuple(np.array(j) + lMPD*(i-1))] = copy.copy(aMPD[j])
# Molecule 2
for i in range(1, nTMC + 1):
    aTMC = copy.copy(TMC.angles)
    for j in aTMC:
        # Add new angle with updated atom indexes
        TOP.angles[tuple(np.array(j) + lTMC*(i-1) + lMPD*nMPD)] = copy.copy(aTMC[j])



### Dihedrals Section ###
# Molecule 1
for i in range(2, nMPD + 1):
    dMPD = copy.copy(MPD.dihedrals)
    for j in dMPD:
        # Add new dihedral with updated atom indexes
        TOP.dihedrals[tuple(np.array(j) + lMPD*(i-1))] = copy.copy(dMPD[j])
# Molecule 2
for i in range(1, nTMC + 1):
    dTMC = copy.copy(TMC.dihedrals)
    for j in dTMC:
        # Add new dihedral with updated atom indexes
        TOP.dihedrals[tuple(np.array(j) + lTMC*(i-1) + lMPD*nMPD)] = copy.copy(dTMC[j])



### Improper Dihedrals Section ###
# Molecule 1
for i in range(2, nMPD + 1):
    iMPD = copy.copy(MPD.improper_dihedrals)
    for j in iMPD:
        # Add new improper dihedral with updated atom indexes
        TOP.improper_dihedrals[tuple(np.array(j) + lMPD*(i-1))] = copy.copy(iMPD[j])
# Molecule 2
for i in range(1, nTMC + 1):
    iTMC = copy.copy(TMC.improper_dihedrals)
    for j in iTMC:
        # Add new improper dihedral with updated atom indexes
        TOP.improper_dihedrals[tuple(np.array(j) + lTMC*(i-1) + lMPD*nMPD)] = copy.copy(iTMC[j])



### Pairs Section ###
# Molecule 1
for i in range(2, nMPD + 1):
    pMPD = copy.copy(MPD.pairs)
    for j in pMPD:
        # Add new pair with updated atom indexes
        TOP.pairs[tuple(np.array(j) + lMPD*(i-1))] = copy.copy(pMPD[j])
# Molecule 2
for i in range(1, nTMC + 1):
    pTMC = copy.copy(TMC.pairs)
    for j in pTMC:
        # Add new pair with updated atom indexes
        TOP.pairs[tuple(np.array(j) + lTMC*(i-1) + lMPD*nMPD)] = copy.copy(pTMC[j])

### Write new TOP ###
TOP.write('polyamide.top')
