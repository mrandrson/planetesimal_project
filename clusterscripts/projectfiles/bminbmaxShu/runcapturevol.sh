#!/bin/bash
#SBATCH -n 1
#SBATCH -o capturevol.out
#SBATCH -e capturevol.err
#SBATCH -p phys-k
python capturevol.py
