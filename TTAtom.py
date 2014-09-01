#!/usr/bin/env python
class Atom:
    def __init__(self,n,x,y,z,atomType,connectivity):
        self.n = n
        self.x = x
        self.y = y
        self.z = z
        self.type = atomType
        self.connectivity = connectivity 
        self.rotate = False