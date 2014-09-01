#!/usr/bin/env python

class Atom:
    """Contains all information for each atom, including atom number (n),
    xyz coordinates (x, y, z), the type of atom in both string and parameter
    file index number (atomType *note this is a list*), the connectivity list
    corresponding to index number within Python (connectivity), and whether or 
    not the atom is to be rotated for the next torsional rotation.
    """
    def __init__(self,n,x,y,z,atomType,connectivity):
        self.n = n
        self.x = x
        self.y = y
        self.z = z
        self.type = atomType
        self.connectivity = connectivity 
        self.rotate = False