#!/bin/bash
#SBATCH -n 16               
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=16 
#SBATCH -o test.out          
#SBATCH -e test.err          
#SBATCH -p phys-k            

parallel_jobs=16             
simulations=2700 

#This is at R_fin=R_sun/10

rm -rf ./outputdata/*.csv

for (( i=1; i<=$simulations; i++ ))
do
  echo "Run #$i"
  v0=$(echo "scale=4; 916.5*($i/$simulations)" | bc -l)
  echo "v0: $v0"

  srun --exclusive -n 1 python3 testparticle.py "$v0" 350000 "$output_file" "$i" &

  if (( $i % $parallel_jobs == 0 )); then
    wait 
  fi
done

wait

