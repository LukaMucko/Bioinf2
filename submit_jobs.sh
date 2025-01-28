#!/bin/bash

k_values=(165 170 175 185 190 195)

# Submit jobs for each k value
for k in "${k_values[@]}"
do
    qsub -v k=$k job.pbs
    echo "Submitted job with k=$k"
done