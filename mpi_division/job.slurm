#!/bin/sh
#BATCH --job-name=ncorps
#SBATCH --time=1:00:00
#SBATCH --partition=medium
## Pour reserver 1 coeur
#SBATCH --nodes=10
#SBATCH --output=sortie.out
#SBATCH --error=error.e

echo $LD_LIBRARY_PATH
/usr/bin/time -o exe_simu.log mpirun --prefix ${NEC_MPI_ROOT} --bind-to-core -np 120 -npernode 12 exe_simu
