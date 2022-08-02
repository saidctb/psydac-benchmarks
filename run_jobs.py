import os

batch_str =\
"""#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.{nprocs}.%j
#SBATCH -e ./tjob.err.{nprocs}.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J job_{nprocs}
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

# Run the program:

"""

filenames       = ['poisson_3d', 'vector_poisson_3d', 'time_harmonic_maxwell_3d']
nnodes          = [1,2,4,8,16,32,64,128]
ntasks_per_node = 32
ncells          = [32,64,128,256]
degrees         = [2,3,4,5]

script_nc_d = 'srun python3 {filename}.py -n {nc} {nc} {nc} -d {d} {d} {d}\n'

os.makedirs('results', exist_ok=True)

for nn in nnodes:
    batch_script = batch_str.format(nprocs=nn*ntasks_per_node,nnodes=nn,ntasks_per_node=ntasks_per_node, mem=150000,time_limit="6:00:00")

    for f in filenames:
        for nc in ncells:
            for d in degrees:
                batch_script += script_nc_d.format(filename=f, nc=nc, d=d)
            
            batch_script += '\n'

    filename = 'job_{nprocs}.sh'.format(nprocs=nn*ntasks_per_node)
    with open(filename,'w') as f:
        f.write(batch_script)

    os.system('sbatch {}'.format(filename))
    

    
    

        
    


