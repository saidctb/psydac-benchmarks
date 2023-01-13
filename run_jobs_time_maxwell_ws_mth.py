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
#SBATCH --cpus-per-task={nthreads2}
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
export OMP_NUM_THREADS={nthreads}
export OMP_PLACES=cores
export OMP_PROC_BIND=close

"""

nnodes          = [7,7*2,7*2**2,7*2**3,7*2**4]
ntasks_per_node = 1
ncells          = [96,128,160]
degrees         = [2,3]
nthreads        = 32

nnodes          = [7*2,7*2**2,7*2**3]
ntasks_per_node = [2,2,2,2]
nthreads        = [16,16,16,16]
ncells0 = [[40,40,20],[40,40,40],[80,40,40],[80,80,40]]
ncells1 = [[44,44,22],[44,44,44],[88,44,44],[88,88,44]]
ncells2 = [[48,48,24],[48,48,48],[96,48,48],[96,96,48]]
ncells = [ncells2]
degrees = [3,4,5]

script_nc_d = 'srun python3 {filename}.py -n {nc0} {nc1} {nc2} -d {d} {d} {d}\n'

f = 'maxwell3d'
for i,(nn,nta,nth) in enumerate(zip(nnodes, ntasks_per_node, nthreads)):
    batch_script = batch_str.format(filename=f, nprocs=nn*nta*nth,nnodes=nn,ntasks_per_node=nta, nthreads=nth,nthreads2=16, mem=150000,time_limit="3:00:00")
    
    for ncj in ncells:
        nc = ncj[i]
        for d in degrees:
            batch_script += script_nc_d.format(filename=f, nc0=nc[0],nc1=nc[1],nc2=nc[2], d=d)

        batch_script += '\n'

    filename = 'job_{filename}_{nprocs}.sh'.format(filename=f, nprocs=nn*nta*nth)
    with open(filename,'w') as file_:
        file_.write(batch_script)

    os.system('sbatch {}'.format(filename))
