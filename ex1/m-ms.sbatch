#!/bin/bash

#SBATCH --job-name=m-MS
#SBATCH --output=m-ms-output.out
#SBATCH --nodes=32
#SBATCH --exclusive
#SBATCH --constraint=ib

module load openmpi
mpiexec -n 2 ./mandelbrot 1000000
mpiexec -n 4 ./mandelbrot 1000000
mpiexec -n 8  ./mandelbrot 1000000
mpiexec -n 16 ./mandelbrot 1000000
mpiexec -n 32 ./mandelbrot 1000000

