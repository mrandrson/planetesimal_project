#!/bin/bash

rm energies.txt
make clean
make

if [ $? -ne 0 ]; then
    echo "Build failed. Exiting."
    exit 1
fi

v0_start=18.33
v0_end=549.9
num_points=2500  

for i in $(seq 0 $num_points); do
    v0=$(awk -v s=$v0_start -v e=$v0_end -v n=$num_points -v i=$i 'BEGIN { print s + i*(e - s)/n }')
    echo "Running simulation with v0 = $v0"
    ./out $v0 350000
done

python energies.py

