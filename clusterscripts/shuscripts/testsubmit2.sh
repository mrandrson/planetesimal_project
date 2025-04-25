#!/bin/bash
#SBATCH -n 16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -o test.out
#SBATCH -e test.err
#SBATCH -p phys-k

parallel_jobs=16
num_v0=100
simulations_per_v0=5
total_simulations=$(( num_v0 * simulations_per_v0 ))

rm -rf ./outputdata/*.csv

index=1

for (( i=1; i<=$num_v0; i++ ))
do
  v0=$(echo "scale=4; 1833*($i/$num_v0)" | bc -l)

  for (( j=1; j<=$simulations_per_v0; j++ ))
  do
    echo "Run #$index (v0: $v0, repetition: $j)"

    srun --exclusive -n 1 python3 testparticle.py "$v0" 350000 "$output_file" "$index" &

    if (( index % parallel_jobs == 0 )); then
      wait
    fi

    ((index++))
  done
done

wait
