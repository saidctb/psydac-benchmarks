#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.h1_projector.27.%j
#SBATCH -e ./tjob.err.h1_projector.27.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J job_h1_projector_27
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=27
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

srun -n 1 python3 h1_projector.py -n 8 8 8 -d 2 2 2
srun -n 1 python3 h1_projector.py -n 8 8 8 -d 3 3 3
srun -n 1 python3 h1_projector.py -n 8 8 8 -d 4 4 4
srun -n 1 python3 h1_projector.py -n 8 8 8 -d 5 5 5


srun -n 8 python3 h1_projector.py -n 16 16 16 -d 2 2 2
srun -n 8 python3 h1_projector.py -n 16 16 16 -d 3 3 3
srun -n 8 python3 h1_projector.py -n 16 16 16 -d 4 4 4
srun -n 8 python3 h1_projector.py -n 16 16 16 -d 5 5 5


srun -n 27 python3 h1_projector.py -n 24 24 24 -d 2 2 2
srun -n 27 python3 h1_projector.py -n 24 24 24 -d 3 3 3
srun -n 27 python3 h1_projector.py -n 24 24 24 -d 4 4 4
srun -n 27 python3 h1_projector.py -n 24 24 24 -d 5 5 5



