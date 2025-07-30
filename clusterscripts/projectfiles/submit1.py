import os

def writeScript(v0, N_exp, t, saveN):
	lines = []
	lines = lines + ['#!/bin/bash\n']
	lines = lines + ['#SBATCH -n 1\n']
	lines = lines + ['#SBATCH -o outputShu.txt\n']
	lines = lines + ['#SBATCH -e errorShu.txt\n']
	lines = lines + ['#SBATCH -p phys-k\n']
	lines = lines + ['python ShuCluster2.py ' + str(v0) + ' ' + str(N_exp) + ' ' + str(t) + ' ' + str(saveN) + '\n']
	arc = open('submit1.sh', 'w')
	arc.writelines(lines)
	arc.close()

def create_output_directories():
	for folder in ['./submit1data', './data/positions']:
		if not os.path.exists(folder):
			os.makedirs(folder)

create_output_directories()


writeScript(20, 10, 350000, './submit1data/s1.txt')
os.system('sbatch submit1.sh')


