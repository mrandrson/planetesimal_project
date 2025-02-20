#!/bin/bash
#SBATCH -n 16               
#SBATCH --nodes=1            
#SBATCH --ntasks-per-node=16 
#SBATCH -o test.out          
#SBATCH -e test.err          
#SBATCH -p phys-k            

parallel_jobs=16             
simulations=1000 

for (( i=1; i<=$simulations; i++ ))
do
  echo "Run #$i"
  v0=$(echo "scale=4; 1833*($i/$simulations)" | bc -l)
  echo "v0: $v0"

  output_file="./outputdata/testpartout$i.h5"

  if [ -f "$output_file" ]; then
    rm -rf "$output_file"
    echo "Deleted file: $output_file"
  fi

  srun --exclusive -n 1 python3 testparticle.py "$v0" 350000 "$output_file" "$i" &

  if (( $i % $parallel_jobs == 0 )); then
    wait 
  fi
done

wait

