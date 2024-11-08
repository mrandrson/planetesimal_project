#!/bin/bash

parallel_jobs=2

for i in {1..250}
do
  echo "Run #$i"
  v0=$(echo "scale=4; 200*($i/250)" | bc -l)
  echo "v0: $v0"
  
  output_file="testpartout$i.h5"
  
  if [ -f "$output_file" ]; then
    rm "$output_file"
    echo "Deleted file: $output_file"
  fi
  
  python3 testparticle.py "$v0" 350000 "$output_file" &

  if (( $i % $parallel_jobs == 0 )); then
    wait
  fi
done

wait

