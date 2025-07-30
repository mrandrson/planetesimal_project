import os
def writeScript(v0, N_exp, t, saveN):
    lines = []
    lines = lines + ['#!/bin/bash\n']
    lines = lines + ['#SBATCH -n 1\n']
    lines = lines + ['#SBATCH -o output.txt\n']
    lines = lines + ['#SBATCH -e error.txt\n']
    lines = lines + ['#SBATCH -p phys-k\n']
    lines = lines + ['python LarsonCluster.py ' + str(v0) + ' ' + str(N_exp) + ' ' + str(t) + ' ' + str(saveN) + '\n']
    arc = open('submit1.sh', 'w')
    arc.writelines(lines)
    arc.close()

for i in range(0, 100):
    writeScript(20*(i+1), 5, -1000, './data1/s' + str(i+1) + '.txt')
    os.system('sbatch submit1.sh')

for i in range(0, 100):
    writeScript(20*(i+1), 5, -10000, './data2/s' + str(i+1) + '.txt')
    os.system('sbatch submit1.sh')

for i in range(0, 100):
    writeScript(20*(i+1), 5, -100000, './data3/s' + str(i+1) + '.txt')
    os.system('sbatch submit1.sh')
