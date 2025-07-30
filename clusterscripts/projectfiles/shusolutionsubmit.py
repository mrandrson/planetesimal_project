import os
import numpy as np

def writeScript(A, Nsteps):
    lines = []
    lines.append('#!/bin/bash\n')
    lines.append('#SBATCH -n 1\n') 
    lines.append(f'#SBATCH -o shusolution.out\n')
    lines.append(f'#SBATCH -e shusolution.err\n')
    lines.append('#SBATCH -p phys-k\n') 
    lines.append('python shuSolution.py {:.12f} {}\n'.format(A, int(Nsteps)))
    script_name = 'shusolution_{:.12f}.sh'.format(A)
    with open(script_name, 'w') as arc: 
        arc.writelines(lines) 
    return script_name 
 
Nsteps = 1e5
 
A_values = 2+np.logspace(-1, -12, 12)

for A in A_values:
    script_name = writeScript(A, Nsteps)
    os.system('sbatch {}'.format(script_name))

