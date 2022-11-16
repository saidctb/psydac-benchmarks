#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.h1_projector.64.%j
#SBATCH -e ./tjob.err.h1_projector.64.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J job_h1_projector_64
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --mem=150000
#
# Wall clock limit:
#SBATCH --time=1:00:00

module load gcc/9
module load openmpi/4
module load anaconda/3/2020.02 mpi4py/3.0.3

# Run the program:

export OMPI_MCA_mpi_warn_on_fork=0
export OMPI_MCA_mpi_yield_when_idle=1

srun python3 h1_projector.py -n 32 32 32 -d 2 2 2
srun python3 h1_projector.py -n 32 32 32 -d 3 3 3
srun python3 h1_projector.py -n 32 32 32 -d 4 4 4
srun python3 h1_projector.py -n 32 32 32 -d 5 5 5

