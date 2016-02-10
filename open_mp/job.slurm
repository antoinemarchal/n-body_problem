#!/bin/sh
#BATCH --job-name=ncorps-omp
#SBATCH --time=1:00:00
#SBATCH --partition=medium
## Pour reserver 1 coeur
#SBATCH --nodes=1
#SBATCH --output=sortie.out

export OMP_NUM_THREADS=12

/usr/bin/time ./exe_omp
