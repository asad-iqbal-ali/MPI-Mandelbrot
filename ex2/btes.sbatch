#!/bin/bash

#SBATCH --job-name=btes
#SBATCH --output=btes-output.out
#SBATCH --nodes=2
#SBATCH --exclusive
#SBATCH --constraint=ib

module load openmpi
mpiexec -n 2 ./btes

