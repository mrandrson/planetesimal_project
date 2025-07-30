#!/bin/bash
#SBATCH -n 80
#SBATCH -o output_master.txt
#SBATCH -e error_master.txt
#SBATCH -p phys-k
module load python
parallel --jobs 80 < all_jobs.txt
