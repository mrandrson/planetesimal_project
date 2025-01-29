#!/bin/bash
#SBATCH -n 16               
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=16 
#SBATCH -o eerr.out          
#SBATCH -e eerr.err          
#SBATCH -p phys-k

parallel_jobs=16
DIRECTORY='/home/rbanderson/projectfiles/shuscripts/outputdata'

i=0
for csvfile in "$DIRECTORY"/*.csv; do
  srun --exclusive -n 1 python3 analyzedat.py "$csvfile" &
  i=$((i + 1))
  
  if (( i % parallel_jobs == 0 )); then
    wait
  fi
done

wait

