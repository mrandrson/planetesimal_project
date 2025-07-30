#!/bin/bash
#SBATCH -n 1
#SBATCH -o test.out
#SBATCH -e test.err
#SBATCH -p phys-k
python testparticle.py 15 93533 350000 ./data3
