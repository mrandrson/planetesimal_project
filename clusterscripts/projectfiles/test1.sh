#!/bin/bash
#SBATCH -n 1
#SBATCH -o test.out
#SBATCH -e test.err
#SBATCH -p phys-k
python testparticle.py 5 31186 350000 ./data1
