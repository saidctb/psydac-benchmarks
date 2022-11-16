import os

batch_str =\
"""#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.{filename}.{nprocs}.%j
#SBATCH -e ./tjob.err.{filename}.{nprocs}.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J job_{filename}_{nprocs}
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes={nnodes}
#SBATCH --ntasks-per-node={ntasks_per_node}
#SBATCH --mem={mem}
#
# Wall clock limit:
#SBATCH --time={time_limit}

module load gcc/9
module load openmpi/4
module load anaconda/3/2020.02 mpi4py/3.0.3
module load h5py-mpi/2.10

# Run the program:

export OMPI_MCA_mpi_warn_on_fork=0
export OMPI_MCA_mpi_yield_when_idle=1

"""

nnodes          = [1,2,2**2,2**3,2**4,2**5,2**6,2**7]
ntasks_per_node = 35
ncells          = [50,75,100,125,150,175,200]
degrees         = [2,3]

script_nc_d = 'srun python3 {filename}.py -n {nc} {nc} {nc} -d {d} {d} {d}\n'

f = 'maxwell3d'
for nn in nnodes:
    batch_script = batch_str.format(filename=f, nprocs=nn*ntasks_per_node,nnodes=nn,ntasks_per_node=ntasks_per_node, mem=150000,time_limit="6:00:00")
    for nc in ncells:
        for d in degrees:
            batch_script += script_nc_d.format(filename=f, nc=nc, d=d)

        batch_script += '\n'

    filename = 'job_{filename}_{nprocs}.sh'.format(filename=f, nprocs=nn*ntasks_per_node)
    with open(filename,'w') as file_:
        file_.write(batch_script)

    os.system('sbatch {}'.format(filename))
