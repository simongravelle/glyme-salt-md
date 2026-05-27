#!/usr/bin/env python
# coding: utf-8

# PEO–TFSI–Cs polymer electrolyte builder in Python for generating bulk molecular dynamics simulation systems, including PEG chain insertion, TFSI anion and Cs⁺ placement, and automatic GROMACS topology (.top) and configuration (.gro) files.
# Author: Simon Gravelle
# License: GNU GPL v3

import numpy as np

import os
os.system("jupyter-nbconvert --to script generatePEOgromacs.ipynb > /dev/null 2>&1")

from peo_generator import PEOgenerator

from utilities import neighborsearch, randomlocation

### desired total number of atoms from the PEO
NAt_des = 6000 

# create PEO
Nsegments = 8
atomsPEG, bondsPEG, anglesPEG, dihedralsPEG, molarmass = PEOgenerator(Nsegments)
atomsPEG[:,3] = atomsPEG[:,3] - np.mean(atomsPEG[:,3])
atomsPEG[:,4] = atomsPEG[:,4] - np.mean(atomsPEG[:,4])
atomsPEG[:,5] = atomsPEG[:,5] - np.mean(atomsPEG[:,5])
Natoms_per_polymer = len(atomsPEG)
nC = np.sum(atomsPEG.T[1] == 1) # number of carbon per Peg
nO = np.sum((atomsPEG.T[1] == 2) | (atomsPEG.T[1] == 4)) # number of oxygen per Peg
nH = np.sum((atomsPEG.T[1] == 3) | (atomsPEG.T[1] == 5)) # number of hydrogen per Peg
Nchains = np.int32(np.floor(NAt_des/Natoms_per_polymer)+1)

# [O]/[XTFSI]=30 (30 oxygène par alcalin)

nTFSI = np.int32((nO * Nchains)/30)

print()
print("Number of segments", Nsegments)
print("Number of chains =", Nchains)
print("Number of atoms =", Nchains*Natoms_per_polymer)
print("Number of TFSI =", nTFSI)
print()

atomsTFSI = np.loadtxt("DATA/TFSI/positions.dat")
bondsTFSI = np.loadtxt("DATA/TFSI/bonds.dat")
anglesTFSI = np.loadtxt("DATA/TFSI/angles.dat")
dihedralsTFSI = np.loadtxt("DATA/TFSI/dihedrals.dat")

cpt_PEO = 0
cpt_TFSI = 0
cpt_Cs = 0
L = np.max([(np.max(atomsPEG[:,3]) - np.min(atomsPEG[:,3])) + 1, 5]) # initial box size - nm
while ((cpt_PEO<Nchains) | (cpt_TFSI<nTFSI) | (cpt_Cs < nTFSI)):

    box = np.array([L, L, L])

    ## initialise matrix
    atoms = np.zeros((100000,7))
    bonds = np.zeros((100000,2))
    angles = np.zeros((100000,3))
    dihedrals = np.zeros((100000,4))
    cpt_atoms = 0
    cpt_bonds = 0
    cpt_angles = 0
    cpt_dihedrals = 0
    cpt_res = 0
    cpt_PEO = 0
    cpt_TFSI = 0
    cpt_Cs = 0
    resName = ["" for x in range(1000000)]
    atoName = ["" for x in range(1000000)]

    fail_attempt = 0
    atoms_PEO = 0
    while (cpt_PEO < Nchains) & (fail_attempt<1e3):

        nO, nH, nC = 0, 0, 0

        insert = 1
        x,y,z = randomlocation(L,L,L)
        if cpt_PEO > 0:
            d = neighborsearch(atoms[:cpt_atoms].T[4:].T, atomsPEG.T[3:].T, x, y, z, box)
            if d < 0.3:
                # Overlap, don't insert
                insert = 0

        if insert == 1:
            for m in bondsPEG:
                bonds[cpt_bonds] = m[0]+cpt_atoms, m[1]+cpt_atoms
                cpt_bonds += 1

            for m in anglesPEG:
                angles[cpt_angles] = m[0]+cpt_atoms, m[1]+cpt_atoms, m[2]+cpt_atoms
                cpt_angles += 1

            for m in dihedralsPEG:
                dihedrals[cpt_dihedrals] = m[0]+cpt_atoms, m[1]+cpt_atoms, m[2]+cpt_atoms, m[3]+cpt_atoms
                cpt_dihedrals += 1

            for m in atomsPEG:
                atomtype = m[1]
                atoms[cpt_atoms] = cpt_atoms+1, cpt_res+1, atomtype, m[2], m[3]+x, m[4]+y, m[5]+z 
                resName[cpt_atoms] = 'PEG'
                if atomtype == 1:
                    nC += 1
                    atoName[cpt_atoms] = 'C'+str(nC)
                elif atomtype == 2:
                    nO += 1
                    atoName[cpt_atoms] = 'O'+str(nO)              
                elif atomtype == 3:
                    nH += 1
                    atoName[cpt_atoms] = 'H'+str(nH)
                elif atomtype == 4:
                    nO += 1
                    atoName[cpt_atoms] = 'O'+str(nO)
                elif atomtype == 5:
                    nH += 1
                    atoName[cpt_atoms] = 'H'+str(nH)
                cpt_atoms += 1
                atoms_PEO += 1
            cpt_res += 1
            cpt_PEO += 1
        else:
            fail_attempt += 1

    atoms_TFSI = 0
    success_attempt = 0
    while (cpt_TFSI < nTFSI) & (fail_attempt<1e3):

        nCi = 0
        nFli = 0
        nSui = 0
        nNii = 0
        nOi = 0

        insert = 1
        x,y,z = randomlocation(L,L,L)
        d = neighborsearch(atoms[:cpt_atoms].T[4:].T, atomsTFSI.T[4:7].T/10, x, y, z, box)
        if d < 0.3: # Overlap, don't insert
            insert = 0

        if insert == 1:
            for m in bondsTFSI[:,2:]:
                bonds[cpt_bonds] = m[0]+cpt_atoms, m[1]+cpt_atoms
                cpt_bonds += 1

            for m in anglesTFSI[:,2:]:
                angles[cpt_angles] = m[0]+cpt_atoms, m[1]+cpt_atoms, m[2]+cpt_atoms
                cpt_angles += 1

            for m in dihedralsTFSI[:,2:]:
                dihedrals[cpt_dihedrals] = m[0]+cpt_atoms, m[1]+cpt_atoms, m[2]+cpt_atoms, m[3]+cpt_atoms
                cpt_dihedrals += 1

            for m in atomsTFSI:
                atomtype = m[2]
                atomcharge = m[3]
                atomtype += 5 # shift for PEG
                atoms[cpt_atoms] = cpt_atoms+1, cpt_res+1, atomtype, atomcharge, m[4]/10+x, m[5]/10+y, m[6]/10+z 
                resName[cpt_atoms] = 'TFSI'
                if atomtype == 6.0:
                    nCi += 1
                    atoName[cpt_atoms] = 'Ci'+str(nCi)
                elif atomtype == 7.0:
                    nFli += 1
                    atoName[cpt_atoms] = 'Fli'+str(nFli)
                elif atomtype == 8.0:
                    nSui += 1
                    atoName[cpt_atoms] = 'Sui'+str(nSui)
                elif atomtype == 9.0:
                    nNii += 1
                    atoName[cpt_atoms] = 'Nii'+str(nNii) 
                elif atomtype == 10.0:
                    nOi += 1
                    atoName[cpt_atoms] = 'Oi'+str(nOi)
                else:
                    print("warning, unknown atoms")
                cpt_atoms += 1
                atoms_TFSI += 1
            cpt_res += 1
            cpt_TFSI += 1
            success_attempt += 1
        else:
            fail_attempt += 1

    fail_attempt = 0
    while (cpt_Cs < nTFSI) & (fail_attempt<1e3):

        insert = 1
        x,y,z = randomlocation(L,L,L)

        d = neighborsearch(atoms[:cpt_atoms].T[4:].T, np.array([[0, 0, 0]]), x, y, z, box)
        if d < 0.3: # Overlap, don't insert
            insert = 0

        if insert == 1:
            cpt_Cs += 1
            atomtype = 10
            atomcharge = 1
            atoms[cpt_atoms] = cpt_atoms+1, cpt_res+1, atomtype, atomcharge, x, y, z 
            resName[cpt_atoms] = 'Cs'
            atoName[cpt_atoms] = 'Cs1'
            cpt_atoms += 1
            cpt_res += 1
        else:
            fail_attempt += 1

    print(str(cpt_PEO)+' PEO')
    print(str(cpt_TFSI)+' TFSI')
    print(str(cpt_Cs)+' Cs')
    print('box size = '+str(np.round(L,1))+' nm')
    L += 1 # nm

atoms = atoms[:cpt_atoms]
bonds = bonds[:cpt_bonds]
angles = angles[:cpt_angles]
dihedrals = dihedrals[:cpt_dihedrals]

assert atoms_PEO == Nchains*Natoms_per_polymer
assert atoms_TFSI == cpt_TFSI * 15
expected_res = Nchains + cpt_Cs + cpt_TFSI
assert expected_res == cpt_res
expected_atom = Nchains*Natoms_per_polymer + cpt_Cs + cpt_TFSI * 15
assert expected_atom == cpt_atoms

f = open('topol.top', 'w')
f.write('#include "ff/forcefield.itp"\n')
f.write('#include "ff/peo.itp"\n')
f.write('#include "ff/tfsi.itp"\n')
f.write('#include "ff/cs.itp"\n\n')
f.write('[ System ]\n')
f.write('bulk PEO-CsTFSI system\n\n')
f.write('[ Molecules ]\n\n')
f.write('PEO '+ str(cpt_PEO)+'\n')
f.write('TFSI '+ str(cpt_TFSI)+'\n')
f.write('Cs '+ str(cpt_Cs)+'\n')
f.close()

f = open('conf.gro', 'w')
f.write('bulk PEO-CsTFSI system\n')
f.write(str(cpt_atoms)+'\n')
nc, no, nh, nm = 0,0,0,0
for n in range(cpt_atoms): 
    f.write("{: >5}".format(np.int32(atoms[n][1]))) # residue number (5 positions, integer) 
    f.write("{: >5}".format(str(resName[n]))) # residue name (5 characters) 
    f.write("{: >5}".format(str(atoName[n]))) # residue name (5 characters) 
    f.write("{: >5}".format(str(np.int32(n+1)))) # atom number (5 positions, integer)
    f.write("{: >8}".format(str("{:.3f}".format(atoms[n][4])))) # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
    f.write("{: >8}".format(str("{:.3f}".format(atoms[n][5])))) # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places) 
    f.write("{: >8}".format(str("{:.3f}".format(atoms[n][6])))) # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places) 
    f.write("\n")
f.write("{: >10}".format(str("{:.5f}".format(L))))
f.write("{: >10}".format(str("{:.5f}".format(L))))
f.write("{: >10}".format(str("{:.5f}".format(L))))
f.write("\n")
f.close()

