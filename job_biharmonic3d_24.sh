#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.ids.out.biharmonic3d.32.%j
#SBATCH -e ./tjob.ids.err.biharmonic3d.32.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J job_biharmonic3d_32
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=150000
#
# Wall clock limit:
#SBATCH --time=3:00:00

module load gcc/9
module load openmpi/4
module load anaconda/3/2020.02 mpi4py/3.0.3
module load h5py-mpi/2.10

# Run the program:

export OMPI_MCA_mpi_warn_on_fork=0
export OMPI_MCA_mpi_yield_when_idle=1

srun -n 1 python3 biharmonic3d_s.py -n 4 4 4 -d 2 2 2
srun -n 1 python3 biharmonic3d_s.py -n 8 8 8 -d 2 2 2
srun -n 8 python3 biharmonic3d_s.py -n 16 16 16 -d 2 2 2
srun python3 biharmonic3d_s.py -n 32 32 32 -d 2 2 2
srun python3 biharmonic3d_s.py -n 32 32 32 -d 2 2 2

srun -n 1 python3 biharmonic3d_s.py -n 4 4 4 -d 3 3 3
srun -n 1 python3 biharmonic3d_s.py -n 8 8 8 -d 3 3 3
srun -n 8 python3 biharmonic3d_s.py -n 16 16 16 -d 3 3 3
srun python3 biharmonic3d_s.py -n 32 32 32 -d 3 3 3
srun python3 biharmonic3d_s.py -n 32 32 32 -d 3 3 3

srun -n 1 python3 biharmonic3d_s.py -n 4 4 4 -d 4 4 4
srun -n 1 python3 biharmonic3d_s.py -n 8 8 8 -d 4 4 4
srun -n 8 python3 biharmonic3d_s.py -n 16 16 16 -d 4 4 4
srun python3 biharmonic3d_s.py -n 32 32 32 -d 4 4 4
srun python3 biharmonic3d_s.py -n 32 32 32 -d 4 4 4
