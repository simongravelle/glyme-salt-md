import numpy as np
from scipy.spatial.transform import Rotation as R
from numpy.linalg import norm
import random

def neighborsearch(neighbor, molecule, x, y, z, box):
    '''Search all neighbor to a molecule in a box and return the closest distance'''
    minr = 10
    for m in molecule:
        x0 = m[0] + x
        y0 = m[1] + y
        z0 = m[2] + z
        dxdydz = np.remainder(neighbor - np.array([x0,y0,z0]) + box/2., box) - box/2.
        minr = np.min([minr,np.min(norm(dxdydz,axis=1))])
    return minr

def randomlocation(Lx,Ly,Lz):
    '''Choose a random location within a given box'''
    txlo, txhi = -Lx/2, Lx/2
    tylo, tyhi = -Ly/2, Ly/2
    tzlo, tzhi = -Lz/2, Lz/2    
    x = random.randint(1,1000)/1000*(txhi-txlo)
    y = random.randint(1,1000)/1000*(tyhi-tylo)
    z = random.randint(1,1000)/1000*(tzhi-tzlo)
    return x, y, z

def randomorientation(XYZ):
    '''3D aleatory rotation of molecule/particule coordinate'''
    rotation_degrees = random.randint(0,9000)/100
    rotation_radians = np.radians(rotation_degrees)
    rotation_axis = np.array([random.randint(0,100)/100, random.randint(0,100)/100, random.randint(0,100)/100])
    rotation_axis /= np.linalg.norm(rotation_axis)
    rotation_vector = rotation_radians * rotation_axis
    rotation = R.from_rotvec(rotation_vector)
    mol_rotated = rotation.apply(XYZ)    
    return mol_rotated.T
