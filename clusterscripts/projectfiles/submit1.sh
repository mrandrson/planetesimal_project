#!/bin/bash
#SBATCH -n 1
#SBATCH -o outputtest.txt
#SBATCH -e errortest.txt
#SBATCH -p phys-k
python testparticle.py 20 10 350 ./test1data/s1.txt
