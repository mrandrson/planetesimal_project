import numpy as np
import os

def writeScript(b_min, b_max, v_0, R_init, b_index, v_index, folder):
    lines = []
    lines = lines + ['#!/bin/bash\n']
    lines = lines + ['#SBATCH -n 1\n']
    lines = lines + ['#SBATCH -o output.txt\n']
    lines = lines + ['#SBATCH -e errorHomogeneous.txt\n']
    lines = lines + ['#SBATCH -p phys-k\n']
    lines = lines + ['python homogeneous.py ' + str(b_min) + ' ' + str(b_max) + ' ' + str(v_0)  + ' ' + str(R_init) + ' ' + str(b_index) + ' ' + str(v_index) + ' ' + str(folder) + '\n']
    arc = open('submit.sh', 'w')
    arc.writelines(lines)
    arc.close()




bvec = np.logspace(9, 17, 80)
vvec = np.logspace(np.log10(18.3), np.log10(18300), 25)
for i in range(0, 79):
    for j in range(0, 25):
         writeScript(bvec[i], bvec[i+1], vvec[j], 1.974*10**15, i, j, 1)
         os.system('sbatch submit.sh')

bvec = np.logspace(9, 17, 80)
for i in range(0, 79):
    for j in range(0, 25):
         writeScript(bvec[i], bvec[i+1], vvec[j], 1.974*10**14, i, j, 2)
         os.system('sbatch submit.sh')

bvec = np.logspace(9, 17, 80)
for i in range(0, 79):
    for j in range(0, 25):
         writeScript(bvec[i], bvec[i+1], vvec[j], 1.974*10**13, i, j, 3)
         os.system('sbatch submit.sh')
