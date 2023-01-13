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

export OMP_NUM_THREADS=1
export OMPI_MCA_mpi_warn_on_fork=0
export OMPI_MCA_mpi_yield_when_idle=1

"""

nnodes          = [1,1,1,1,1,1,2,4,8,16,32,64]
ntasks_per_node = [1,2,4,8,16,32,32,32,32,32,32,32]
ncells0 = [[8,8,8],[16,10,10],[16,16,10],[16,16,16],[32,16,16],[32,32,16],[32,32,32],[64,32,32],[64,64,32],[64,64,64],[128,64,64],[128,128,64]]
ncells1 = [[10,10,10],[20,10,10],[20,20,10],[20,20,20],[40,20,20],[40,40,20],[40,40,40],[80,40,40],[80,80,80],[160,80,80],[160,160,80]]
ncells = [ncells0,ncells1]
degrees = [2,3,4]

script_nc_d = 'srun python3 {filename}.py -n {nc0} {nc1} {nc2} -d {d} {d} {d}\n'

f = 'maxwell3d_ws'
for i,(nn,nt) in enumerate(zip(nnodes, ntasks_per_node)):
    batch_script = batch_str.format(filename=f, nprocs=nn*nt,nnodes=nn,ntasks_per_node=nt, mem=150000,time_limit="5:00:00")
    for ncj in ncells:
        nc = ncj[i]
        for d in degrees:
            batch_script += script_nc_d.format(filename=f, nc0=nc[0],nc1=nc[1],nc2=nc[2], d=d)

        batch_script += '\n'

    filename = 'job_{filename}_{nprocs}.sh'.format(filename=f, nprocs=nn*nt)
    with open(filename,'w') as file_:
        file_.write(batch_script)

    os.system('sbatch {}'.format(filename))
