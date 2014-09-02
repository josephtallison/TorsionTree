#!/usr/bin/env python
import numpy as np
from numpy import linalg as la
import itertools as it
import subprocess
from multiprocessing import Pool
import TTFiles
import TTAtom
import TTTor
import TTEnergy
from TTSysInfo import SysInfo
        
def check_redundancy(bondList,bond):
    for checkBond in bondList:
        if checkBond[0] in bond and checkBond[1] in bond:
            return True
    return False
        
def assign_torsions(atoms):
    """First define the non-hydrogen participating bonds then check for any 
    non-hydrogen atoms connected to them.
    """
    #Define central bond first
    torsions = []
    ccBonds = []
    for atom in atoms:
        if "H" not in atom.type[0]:
            for connectedAtomIndex in atom.connectivity:
                connectedAtom = atoms[connectedAtomIndex]
                if "H" not in connectedAtom.type[0]:
                    if not check_redundancy(ccBonds,[atom,connectedAtom]):
                        ccBonds.append([atom,connectedAtom])
    #Check and see if an additional non H is attached to each atom in the bond                    
    for bond in ccBonds:
        connectedAtoms = []
        for atom in bond:
            for connectedAtomIndex in atom.connectivity:
                connectedAtom = atoms[connectedAtomIndex]
                if "H" not in connectedAtom.type[0]:
                    if connectedAtom not in bond:
                        connectedAtoms.append(connectedAtom)
                        break
        if len(connectedAtoms) == 2:
            torAtoms = [connectedAtoms[0],bond[0],bond[1],connectedAtoms[1]]                  
            torsions.append(TTTor.Torsion(torAtoms))
    return torsions

def check_deltaPhi(deltaPhi):
    """Makes sure that deltaPhi can divide into 360 evenly.  Returns BOOL."""
    status = True
    if (360%deltaPhi != 0):
        status = False
    return status
    
def check_deltaPhi_amount(sysInfo):
    """Makes sure that the amount of torsions that would result from the angle
    step size and the size of the molecule is not above some intractable
    threshold
    """
    amount = (360/sysInfo.deltaPhi)**sysInfo.nTors
    if amount > 10000:
        return False,amount
    else:
        return True,amount
    
def assign_rotating_atoms(atom1,atom2,atoms):
    """Defines which atoms need to be rotated (atom1-->atom2 defines which 
    direction in the chain to assign rotatable atoms).  Sets the atom.rotate 
    logical for use in rotation.
    """      
    atomsToRotate = [atom2]
    for atom in atomsToRotate:
        atom.rotate = True
        for connectedAtomIndex in atom.connectivity:
            connectedAtom = atoms[connectedAtomIndex]
            if ((connectedAtom not in atomsToRotate) and
                    (connectedAtom != atom1)):
                connectedAtom.rotate = True
                atomsToRotate.append(connectedAtom)

def prepare_rotation_matrix(rotationVector):
    """Generates YZ rotation matrix from rotationVector."""
    rotationVectorLength = la.norm(rotationVector)
    cost=rotationVector[2]/rotationVectorLength #direction cos z
    sint=np.sqrt(np.abs(1.0-cost**2))
    if sint == 0:
        cosp=1.0
        sinp=0.0
    else:
        cosp=rotationVector[0]/(rotationVectorLength*sint) #spherical coordinate identities
        sinp=rotationVector[1]/(rotationVectorLength*sint)
    #Y & Z rotation matrices multiplied together
    YZRotationMat=np.matrix([[cost*cosp,cost*sinp,-sint],
                            [-sinp,cosp,0.0],
                            [sint*cosp,sint*sinp,cost]])
    return YZRotationMat

def rotate_atoms(torsion,atoms):
    """Rotates all atoms corresponding to the rotation of torsion by rotating 
    the selected atoms into the YZ plane, performing the rotation about the Z 
    axis, then rotating back to its original orientation.  Sets the x, y, & z
    coordinates within the atom object.
    """
    #Set up easily read variables
    atom1Coord = np.array([torsion.tor_atoms[1].x,
                           torsion.tor_atoms[1].y,
                           torsion.tor_atoms[1].z])
    atom2Coord = np.array([torsion.tor_atoms[2].x,
                           torsion.tor_atoms[2].y,
                           torsion.tor_atoms[2].z])
    rotPhi = torsion.newAngle + torsion.currentAngle()
    rotationVector = atom2Coord - atom1Coord

    #Create Final Rotation Matrix
    YZRotMat = prepare_rotation_matrix(rotationVector)
    ZRotMat=np.matrix([[np.cos(rotPhi),-np.sin(rotPhi),0.0],
                      [np.sin(rotPhi),np.cos(rotPhi),0.0],
                      [0.0,0.0,1.0]])
    rotationMatrix = np.dot(YZRotMat.T,np.dot(ZRotMat,YZRotMat))
        
    #Perform Rotation
    for atom in atoms:
        if atom.rotate:
            atomToRotateCoord = np.array([atom.x,atom.y,atom.z])
            rotationVector = np.matrix(atomToRotateCoord - atom1Coord)
            atomToRotateCoord = np.dot(rotationVector,rotationMatrix.T)           
            atomToRotateCoord = atomToRotateCoord + atom1Coord
            atom.x = float(atomToRotateCoord[[0],[0]])
            atom.y = float(atomToRotateCoord[[0],[1]])
            atom.z = float(atomToRotateCoord[[0],[2]])
        
def check_close_contacts(atoms):
    """Checks distance between atoms pairwise.  Of the pairwise matrix, only
    the upper diagonal is evaluated.  Returns True if the atoms are too close,
    otherwise it will return False.
    """
    for i in xrange(len(atoms)):
        for j in xrange(i+1,len(atoms)):
            xDif = (atoms[i].x-atoms[j].x)**2
            yDif = (atoms[i].y-atoms[j].y)**2
            zDif = (atoms[i].z-atoms[j].z)**2
            totDif = np.sqrt(xDif + yDif + zDif)
            if totDif < 0.5:    
                return True
    return False
 
def create_possible_torsions(deltaPhi):
    """Given an angle step-size, returns a list of all possible torsions from 
    180 to -180 (not including -180).
    """
    amount = int(360/deltaPhi)
    possibleValues = []
    possibleValue = 180
    for i in xrange(amount):
        possibleValues.append(possibleValue)
        possibleValue = possibleValue - deltaPhi    
    return possibleValues

def generate_torsion_combinations(possibleTorsions,size):
    for combination in it.product(possibleTorsions,repeat=size):
        yield TTTor.TorsionCombination(combination) 
    
def main_wrapper(torComb):
    #Set angle in torsions object
    for i in xrange(len(torsions)):
        torsions[i].newAngle = np.radians(torComb.combination[i])
        
    #Check which torsions need to be rotated
    for torsion in torsions:
        if torsion.newAngle != torsion.currentAngle():
        
            #Set atoms.rotate to false
            for atom in atoms:
                atom.rotate = False
            
            #Check which atoms need to be rotated & setting atom.rotate
            assign_rotating_atoms(torsion.tor_atoms[1],
                                  torsion.tor_atoms[2],
                                  atoms)
            #Rotate atoms
            rotate_atoms(torsion,atoms)
                
    #Check if atoms are too close together
    if check_close_contacts(atoms):
        return
        
    #Write file for energy evaluation
    torComb.filename = "new_octane"
    name = TTFiles.make_filename(torComb)
    TTFiles.write_xyz(atoms,torComb.filename)
    
    #Perform single point energy evaluation
    torComb.energy = TTEnergy.get_single_point_energy(torComb.filename)
    
    #Check to see if single point energy is above the threshold for minimization
    if torComb.energy > sysInfo.max_energy():
        subprocess.Popen(["rm",torComb.filename])
        return
        
    #Perform minimization
    torComb.min_energy = TTEnergy.get_minimum_energy(torComb.filename)
    print torComb.filename,torComb.min_energy
    
if __name__ == '__main__':
    #SET UP DATASTRUCTURES    
    #Generate an array of atom objects
    sysInfo = SysInfo(raw_input("Filename of Tinker XYZ file: "))
    atoms = TTFiles.get_xyz(sysInfo.startFile)
    sysInfo.startEnergy = TTEnergy.get_single_point_energy(sysInfo.startFile)

    #Generate an array of torsion objects
    torsions = assign_torsions(atoms)    

    #Populate system info & list of possible torsion angles
    sysInfo.nTors = len(torsions)
    
    while (not sysInfo.deltaPhiStatus):
        sysInfo.deltaPhi = float(raw_input("Delta Phi - Angle step size:  "))
        sysInfo.deltaPhiStatus = check_deltaPhi(sysInfo.deltaPhi)
        if not sysInfo.deltaPhiStatus:
            print "\n------------------------------------------"
            print "Delta phi does not divide evenly into 360."
            print "Please provide another value for Delta Phi"
            print "------------------------------------------"
        sysInfo.deltaPhiStatus, sysInfo.numberOfConformers = \
                                check_deltaPhi_amount(sysInfo)
        if not sysInfo.deltaPhiStatus:
            print "\n--------------------------------------------------------"
            print ("The current delta phi will result in %s" % 
                   sysInfo.numberOfConformers)
            print "conformations.  Please provide a larger angle step size." 
            print "--------------------------------------------------------"     
    sysInfo.possibleTorsions = create_possible_torsions(sysInfo.deltaPhi)
    sysInfo.pool = Pool(6)

    #PERFORM SEARCH
    #Iterate over every possible combination of torsion angles & minimize
    sysInfo.pool.map(main_wrapper,
                     generate_torsion_combinations(sysInfo.possibleTorsions,
                                                    sysInfo.nTors))





















