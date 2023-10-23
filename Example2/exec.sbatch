#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -J IC2
#SBATCH -t 240:00:00
#SBATCH -o IC2.out
#SBATCH -e IC2.err

ml apps/matlab
matlab -r "SbatchFunction()"
