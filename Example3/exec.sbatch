#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -J shalWa_500
#SBATCH -t 24:00:00
#SBATCH -o out_500.out
#SBATCH -e err_500.err

python3 ex3_artificialDynamics.py;