#!/bin/bash

k_values=(1100, 1200, 1300, 1400, 1600, 1700, 1800, 1900)

# Submit jobs for each k value
for k in "${k_values[@]}"
do
    qsub -v k=$k job.pbs
    echo "Submitted job with k=$k"
done