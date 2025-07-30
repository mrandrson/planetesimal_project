#!/bin/bash

for job_start in 0 1000 2000; do
    job_end=$((job_start + 999))
    if [ $job_end -gt 2399 ]; then
        job_end=2399
    fi
    sbatch --array=${job_start}-${job_end}%100 submit.sh
done

