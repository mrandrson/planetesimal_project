#!/bin/bash
#SBATCH -n 1
#SBATCH -o capturevol.out          
#SBATCH -e capturevol.err          
#SBATCH -p phys-k 

for i in {1..10}
do
	python3 getCaptureVolume.py $i
done
