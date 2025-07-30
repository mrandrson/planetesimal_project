import numpy as np
import os
import math

'''
def writeScript(bmin, bmax, v_0, t_final, bindex, vindex, folder, reps):
    lines = []
    lines = lines + ['#!/bin/bash\n']
    lines = lines + ['#SBATCH -n 1\n']
    lines = lines + ['#SBATCH -o output.txt\n']
    lines = lines + ['#SBATCH -e errorShu.txt\n']
    lines = lines + ['#SBATCH -p phys-k\n']
    lines = lines + ['python Shu.py ' + str(bmin) + ' ' + str(bmax) + ' '  + str(v_0) + ' ' + str(t_final) +  ' '  + str(bindex) + ' ' + str(vindex) + ' ' + str(folder) + ' ' + str(reps) + '\n']
    arc = open('submit.sh', 'w')
    arc.writelines(lines)
    arc.close()

cs = math.sqrt(3.36)*100
bvec = np.logspace(9, 17, 80)
Nv = 30
vvec = np.logspace(np.log10(cs/10), np.log10(cs*5), Nv)
for i in range(0, 30):
    for j in range(0, Nv):
         writeScript(bvec[i], bvec[i+1],  vvec[j], 1.0829*10.0**13, i, j, 1, 100)
         os.system('sbatch submit.sh')
'''
'''
bvec = np.logspace(9, 17, 80)
vvec = np.logspace(np.log10(cs/10), np.log10(cs*10), Nv)
for i in range(0, 79):
    for j in range(0, Nv):
         writeScript(bvec[i], bvec[i+1],  vvec[j], 1.0829*10**12, i, j, 2, 100)
         os.system('sbatch submit.sh')

bvec = np.logspace(9, 17, 80)
vvec = np.logspace(np.log10(cs/10), np.log10(cs*10), Nv)
for i in range(0, 79):
    for j in range(0, Nv):
         writeScript(bvec[i], bvec[i+1],  vvec[j], 1.0829*10.0**13, i, j, 3, 100)
         os.system('sbatch submit.sh')
'''
'''
cs = math.sqrt(3.36)*100
bvec = np.logspace(9, 17, 80)
Nv = 30
vvec = np.logspace(np.log10(cs/10), np.log10(cs*5), Nv)
for i in range(0, 79):
    for j in range(0, Nv):
         writeScript(bvec[i], bvec[i+1],  vvec[j], 1.0829*10.0**10, i, j, 4, 100)
         os.system('sbatch submit.sh')


cs = math.sqrt(3.36)*100
bvec = np.logspace(9, 17, 80)
Nv = 30
vvec = np.logspace(np.log10(cs/10), np.log10(cs*5), Nv)
for i in range(0, 79):
    for j in range(0, Nv):
         writeScript(bvec[i], bvec[i+1],  vvec[j], 1.0829*10.0**9, i, j, 5, 100)
         os.system('sbatch submit.sh')


cs = math.sqrt(3.36)*100
bvec = np.logspace(9, 17, 80)
Nv = 30
vvec = np.logspace(np.log10(cs/10), np.log10(cs*5), Nv)
for i in range(0, 79):
    for j in range(0, Nv):
         writeScript(bvec[i], bvec[i+1],  vvec[j], 1.0829*10.0**8, i, j, 6, 100)
         os.system('sbatch submit.sh')
'''

import os
import numpy as np
import math

def writeMasterScript(bvec, vvec, t_final, folder, reps, max_jobs=30):
    script_name = 'submit_master.sh'
    job_list_name = 'all_jobs.txt'

    with open(job_list_name, 'w') as job_list:
        for i in range(len(bvec) - 1):
            for j, v_0 in enumerate(vvec):
                job_list.write(f'python ShuOutput.py {bvec[i]} {bvec[i+1]} {v_0} {t_final} {i} {j} {folder} {reps}\n')

    with open(script_name, 'w') as arc:
        arc.writelines([
            '#!/bin/bash\n',
            f'#SBATCH -n {max_jobs}\n',
            '#SBATCH -o output_master.txt\n',
            '#SBATCH -e error_master.txt\n',
            '#SBATCH -p phys-k\n',
            'module load python\n',
            f'parallel --jobs {max_jobs} < {job_list_name}\n'
        ])

    return script_name

cs = math.sqrt(3.36) * 100
bvec = np.logspace(9, 17, 80)
Nv = 30
vvec = np.logspace(np.log10(cs / 10), np.log10(cs * 5), Nv)

script = writeMasterScript(bvec, vvec, 1.0829 * 10.0**13, 1, 100, max_jobs=80)
os.system(f'sbatch {script}')

