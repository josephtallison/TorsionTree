#!/usr/bin/env python
import numpy as np
from numpy import linalg as la

class Torsion:
    def __init__(self,tor_atoms):
        self.tor_atoms = tor_atoms
        self.newAngle = None
    def currentAngle(self):
        atom0 = np.array([self.tor_atoms[0].x,
                          self.tor_atoms[0].y,
                          self.tor_atoms[0].z])
        atom1 = np.array([self.tor_atoms[1].x,
                          self.tor_atoms[1].y,
                          self.tor_atoms[1].z])
        atom2 = np.array([self.tor_atoms[2].x,
                          self.tor_atoms[2].y,
                          self.tor_atoms[2].z])
        atom3 = np.array([self.tor_atoms[3].x,
                          self.tor_atoms[3].y,
                          self.tor_atoms[3].z])
        v1 = atom1 - atom0
        v2 = atom2 - atom1
        v3 = atom3 - atom2
        v1v2 = np.cross(v2,v1)
        v2v3 = np.cross(v3,v2)
        v1v2_m = la.norm(v1v2)
        v2v3_m = la.norm(v2v3)
        angle = np.dot(v1v2,v2v3)/(v1v2_m*v2v3_m)
        if angle < -1.0:
            angle = -1.0
        if angle > 1.0:
            angle = 1.0
        angle = np.arccos(angle)
        test = np.cross(v1v2,v2v3)
        test = np.dot(test,v2)
        if test > 0:
            angle = - angle    
        return angle
        
class TorsionCombination:
    def __init__(self,combination):
        self.combination = combination
        self.minCombination = None
        self.energy = None