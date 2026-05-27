#!/usr/bin/env python
# coding: utf-8

# PEO polymer generator in Python for constructing molecular structures and generating GROMACS topology files (.itp), including atoms, bonds, angles, and dihedrals for molecular dynamics simulations.
# Author: Simon Gravelle
# License: GNU GPL v3

import numpy as np
import warnings
warnings.filterwarnings('ignore')

def PEOgenerator(Nseg, Lx=40, Ly=40, Lz=40):
    endpatch = np.loadtxt('DATA/endpatch_CH3.dat')
    monomer = np.loadtxt('DATA/monomer.dat')
    v = 0.14 # distance between 2 monomers

    atoms = []
    cptatoms = 0
    shift = 0
    for id, type, q, x, y, z in endpatch:
        atoms.append([cptatoms+1, type, q, -x, y, z])
        cptatoms+= 1
    shift = np.max(np.array(atoms)[:,3]) + v
    # place N monomers
    for seg in range(Nseg-1):
        for id, type, q, x, y, z in monomer:
            atoms.append([cptatoms+1, type, q, x+shift-np.min(monomer[:,3]), y, z])
            cptatoms += 1  
        shift = np.max(np.array(atoms)[:,3]) + v
    for id, type, q, x, y, z in endpatch:
        atoms.append([cptatoms+1, type, q, x+shift+np.min(monomer[:,3]), y, z])
        cptatoms+= 1
    shift = np.max(np.array(atoms)[:,3]) + v

    atoms = np.array(atoms)
    car = atoms[atoms.T[1] == 1]
    hyd = atoms[(atoms.T[1] == 3) | (atoms.T[1] == 5)]
    oxy = atoms[(atoms.T[1] == 2) | (atoms.T[1] == 4)]

    atoms.T[3] -= np.mean(atoms.T[3])
    atoms.T[3] += Lx/2
    atoms.T[4] -= np.mean(atoms.T[4])
    atoms.T[4] += Ly/2
    atoms.T[5] -= np.mean(atoms.T[5])
    atoms.T[5] += Lz/2

    mC = 12 # m/mol
    mO = 16 # g/mol
    mH = 1 # g/mol
    molmass = len(car)*mC+len(oxy)*mO+len(hyd)*mH
    print('PEO - '+str(molmass)+' g/mol')

    bonds = []
    cptbonds = 0
    # carbon - carbon bonds between monomers
    cpt_CC = 0
    for C1 in car:
        id1, _, _, x1, y1, z1 = C1
        for C2 in car:
            id2, _, _, x2, y2, z2 = C2
            if id1 < id2:
                d = np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
                if d < 0.25:
                    bonds.append([id1, id2])
                    cptbonds += 1
                    cpt_CC += 1
    assert cpt_CC == (Nseg-2)+2

    cpt_CO = 0
    # carbon - oxygen bonds
    xyz = car.T[3:].T
    for n0 in range(len(oxy)):
        xyz0 = oxy[n0][3:]
        idx0 = np.int32(oxy[n0][0])
        d = np.sqrt((xyz.T[0]-xyz0[0])**2+(xyz.T[1]-xyz0[1])**2+(xyz.T[2]-xyz0[2])**2)
        where = np.where((d > 0) & (d < 0.15))
        for w in where[0]:
            idx1 = np.int32(car[w][0])
            if idx0<idx1:
                bonds.append([idx0, idx1])
            else:
                bonds.append([idx1, idx0])
            cptbonds += 1
            cpt_CO += 1
    assert cpt_CO == len(oxy)*2

    # carbon - hydrogen bonds
    cpt_CH = 0
    xyz = car.T[3:].T
    for n0 in range(len(hyd)):
        xyz0 = hyd[n0][3:]
        idx0 = np.int32(hyd[n0][0])
        d = np.sqrt((xyz.T[0]-xyz0[0])**2+(xyz.T[1]-xyz0[1])**2+(xyz.T[2]-xyz0[2])**2)
        where = np.where((d > 0) & (d < 0.11))[0]
        if where.shape == (1,):
            idx1 = car[where][0][0]
            if idx0<idx1:
                bonds.append([idx0, idx1])
            else:
                bonds.append([idx1, idx0])
            cptbonds += 1  
            cpt_CH += 1
    assert cpt_CH == len(hyd)

    # oxygen - hydrogen bonds
    xyz = oxy.T[3:].T
    for n0 in range(len(hyd)):
        xyz0 = hyd[n0][3:]
        idx0 = np.int32(hyd[n0][0])
        d = np.sqrt((xyz.T[0]-xyz0[0])**2+(xyz.T[1]-xyz0[1])**2+(xyz.T[2]-xyz0[2])**2)
        where = np.where((d > 0) & (d < 0.11))[0]
        if where.shape == (1,):
            idx1 = oxy[where][0][0]
            if idx0<idx1:
                bonds.append([idx0, idx1])
            else:
                bonds.append([idx1, idx0])
            cptbonds += 1  
            print("Irrelevant for PEO") 
    # remove excess lines and reorder
    bonds = np.array(bonds)
    bonds = bonds[bonds[:, 0].argsort()]    

    angles = np.zeros((10000,3))
    cptangles = 0
    bonded_a = np.append(bonds.T[0],bonds.T[1])
    for a in atoms:
        ida = np.int32(a[0])
        tpa = np.int32(atoms[atoms.T[0] == ida].T[1])[0]
        occurence = np.sum(bonded_a == ida)
        if occurence > 1: # the atom has 2 or more atoms
            id_neighbors = np.unique(bonds[(bonds.T[0] == ida) | (bonds.T[1] == ida)].T[:2].T)
            for idb in id_neighbors:
                for idc in id_neighbors:
                    if (idb != ida) & (idc != ida) & (idb < idc): # avoid counting same angle twice
                        angles[cptangles] = idb, ida, idc
                        cptangles += 1       
    angles = angles[:cptangles]

    dihedrals = np.zeros((10000,4))
    cptdihedrals = 0
    central_angled_a = angles.T[1]
    edge_angled_a = np.append(angles.T[0],angles.T[2])
    for a in atoms:
        ida = np.int32(a[0])
        tpa = np.int32(atoms[atoms.T[0] == ida].T[1])[0]
        if (tpa == 1) | (tpa == 2) | (tpa == 4): # ignore hydrogen
            id_first_neighbor = np.unique(angles[(angles.T[1] == ida)].T[:3].T)
            id_first_neighbor = id_first_neighbor[id_first_neighbor != ida]
            for idb in id_first_neighbor:
                id_second_neighbor = np.unique(angles[(angles.T[1] == idb)].T[:3].T)
                if len(id_second_neighbor)>0:
                    id_second_neighbor = id_second_neighbor[id_second_neighbor != idb]
                    id_second_neighbor = id_second_neighbor[id_second_neighbor != ida]
                    for idc in id_first_neighbor:
                        if idc != idb:
                            for ide in id_second_neighbor:
                                tpc = np.int32(atoms[atoms.T[0] == idc].T[1])[0]
                                tpe = np.int32(atoms[atoms.T[0] == ide].T[1])[0]
                                if (ida < idb) & (tpc != 3) & (tpe != 3) : 
                                    dihedrals[cptdihedrals] = idc, ida, idb, ide
                                    cptdihedrals += 1
    dihedrals = dihedrals[:cptdihedrals]

    f = open('ff/peo.itp', 'w')
    f.write('[ moleculetype ]\n')
    f.write('PEO   2\n\n')
    f.write('[ atoms ]\n')
    nc = 0
    no = 0
    nh = 0
    for n in range(cptatoms):
        f.write("{: >5}".format(str(n+1))) # atom number
        if atoms.T[1][n] == 1:
            f.write("{: >8}".format('CC32A'))
        elif atoms.T[1][n] == 2:
            f.write("{: >8}".format('OC30A'))
        elif atoms.T[1][n] == 3:
            f.write("{: >8}".format('HCA2'))
        elif atoms.T[1][n] == 4:
            f.write("{: >8}".format('OC311'))
        elif atoms.T[1][n] == 5:
            f.write("{: >8}".format('HCP1'))
        else:
            print('extra atoms')    
        f.write("{: >8}".format(str(1))) # residue number
        f.write("{: >8}".format('PEO')) # residue number
        if atoms.T[1][n] == 1:
            nc += 1
            f.write("{: >8}".format('C'+str(nc))) # atom name
        elif (atoms.T[1][n] == 3) | (atoms.T[1][n] == 5):
            nh += 1
            f.write("{: >8}".format('H'+str(nh))) # atom name
        elif (atoms.T[1][n] == 2) | (atoms.T[1][n] == 4):
            no += 1
            f.write("{: >8}".format('O'+str(no))) # atom name
        f.write("{: >8}".format(str(np.int32(n+1))))
        f.write("{: >8}".format(str("{:.3f}".format(atoms.T[2][n]))))
        if atoms.T[1][n] == 1:
            f.write("{: >8}".format(str("{:.3f}".format(12.011))))
        elif (atoms.T[1][n] == 3) | (atoms.T[1][n] == 5):
            f.write("{: >8}".format(str("{:.3f}".format(1.008))))    
        elif (atoms.T[1][n] == 2) | (atoms.T[1][n] == 4):
            f.write("{: >8}".format(str("{:.3f}".format(15.9994)))) 
        f.write("\n") 
    f.write("\n")  
    f.write('[ bonds ]\n')  
    for n in range(cptbonds):
        f.write("{: >5}".format(str(np.int32(bonds[n][0]))))
        f.write("{: >5}".format(str(np.int32(bonds[n][1]))))
        f.write("{: >5}".format(str(np.int32(1))))
        f.write("\n")
    f.write("\n")  
    f.write('[ angles ]\n')  
    for n in range(cptangles):
        f.write("{: >5}".format(str(np.int32(angles[n][0]))))
        f.write("{: >5}".format(str(np.int32(angles[n][1]))))
        f.write("{: >5}".format(str(np.int32(angles[n][2]))))
        f.write("{: >5}".format(str(np.int32(5))))
        f.write("\n")
    f.write("\n")  
    f.write('[ dihedrals ]\n')  
    for n in range(cptdihedrals):
        f.write("{: >5}".format(str(np.int32(dihedrals[n][0]))))
        f.write("{: >5}".format(str(np.int32(dihedrals[n][1]))))
        f.write("{: >5}".format(str(np.int32(dihedrals[n][2]))))
        f.write("{: >5}".format(str(np.int32(dihedrals[n][3]))))
        f.write("{: >5}".format(str(np.int32(9))))
        f.write("\n")
    f.close()

    return atoms, bonds, angles, dihedrals, molmass

