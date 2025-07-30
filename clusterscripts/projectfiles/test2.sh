#!/bin/bash
#SBATCH -n 1
#SBATCH -o test.out
#SBATCH -e test.err
#SBATCH -p phys-k
python testparticle.py 10 62358 350000 ./data2
