#!/usr/bin/env python
import TTFiles
import subprocess

def get_single_point_energy(atoms,torComb):
    torComb.filename = "new_octane"
    name = TTFiles.make_filename(torComb)
    TTFiles.write_xyz(atoms,torComb.filename)
    command = ("$TINKER/analyze "+torComb.filename+
               " $TINKER/../params/mm3.prm E")
    output = subprocess.check_output(command,shell=True)
    output = output.split()
    for index, out in enumerate(output):
        if "Potential" in out:
            torComb.energy = float(output[index+3])
            