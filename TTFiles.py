#!/usr/bin/env python
from TTAtom import Atom
import csv

def get_xyz(filename):
    """Given an Tinker XYZ coordinate file, this function will extract all of
    the information and place it all within an Atom object.
    """
    atoms = []
    count = 0
    for line in csv.reader(open(filename), delimiter=" ",
                           skipinitialspace=True):
        if count > 0:
            atoms.append(Atom(
                              count,
                              float(line[2]),
                              float(line[3]),
                              float(line[4]),
                              [line[1],int(line[5])],
                              [a-1 for a in map(int,line[6:])]
                              )
                         )
        count += 1
    return atoms

def write_xyz(atoms,filename):
    """Will write a Tinker XYZ coordinate file named "filename" based on the 
    information stored within each Atom class within the atoms array
    """
    outfile = open(filename,"w")
    outfile.write("%6s" % len(atoms)+"\n")
    for atom in atoms:
        writeLine = "%6s" % str(atom.n)
        writeLine = writeLine + "%3s" % atom.type[0]
        writeLine = writeLine + "%14s" % str("{:.6f}".format(atom.x))
        writeLine = writeLine + "%12s" % str("{:.6f}".format(atom.y))
        writeLine = writeLine + "%12s" % str("{:.6f}".format(atom.z))
        writeLine = writeLine + "%6s" % str(atom.type[1])
        for connectedAtomIndex in atom.connectivity:
            writeLine = writeLine + "%6s" % str(connectedAtomIndex+1)
        outfile.write(writeLine+"\n")
    outfile.close()
        
def make_filename(torComb):
    for i in torComb.combination:
        if i < 0:
            lead = "n"
        else:
            lead = ""
        frag = lead + str("%0.0f" % abs(i))
        torComb.filename = torComb.filename + "_" + frag
    torComb.filename = torComb.filename + ".xyz"