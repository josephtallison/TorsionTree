#!/usr/bin/env python

class SysInfo:
    """Contains miscellaneous information regarding the system"""
    def __init__(self,filename):
        self.startFile = filename
        self.startEnergy = None
        self.nTors = None
        self.possibleTors = None
        self.deltaPhi = None
        self.deltaPhiStatus = False
        self.pool = None
        self.numberOfConformers = None
    def max_energy(self):
        max_energy = self.startEnergy+100
        return max_energy
