#!/bin/bash

#SBATCH --job-name=mandelbrot32
#SBATCH --output=m-def-output.out
#SBATCH --nodes=32
#SBATCH --exclusive
#SBATCH --constraint=ib

module load openmpi
mpiexec -n 1 ./mandelbrot
mpiexec -n 2 ./mandelbrot
mpiexec -n 4 ./mandelbrot
mpiexec -n 8  ./mandelbrot
mpiexec -n 16 ./mandelbrot
mpiexec -n 32 ./mandelbrot

