#!/bin/bash
#SBATCH -n 40               
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=40 
#SBATCH -o test.out          
#SBATCH -e test.err          
#SBATCH -p phys-k            

parallel_jobs=40             
simulations=1 

#This is at R_fin=R_sun

rm -rf ./outputdata/simdata.h5

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

