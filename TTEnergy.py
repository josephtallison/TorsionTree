#!/usr/bin/env python
import TTFiles
import subprocess

def get_single_point_energy(filename):
    command = ("$TINKER/analyze "+filename+
               " $TINKER/../params/mm3.prm E")
    output = subprocess.check_output(command,shell=True)
    output = output.split()
    for index, out in enumerate(output):
        if "Potential" in out:
            energy = float(output[index+3])
    return energy
            
def get_minimum_energy(filename):
    command = ("$TINKER/minimize "+filename+
               " $TINKER/../params/mm3.prm 0.01")
    output = subprocess.check_output(command,shell=True)
    output = output.split()
    for index, out in enumerate(output):
        if "Function" in out:
            energy = float(output[index+3])
    return energy
            