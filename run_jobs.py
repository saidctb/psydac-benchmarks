import os

batch_str =\
"""#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.{filename}.{mapping}.{nprocs}.%j
#SBATCH -e ./tjob.err.{filename}.{mapping}.{nprocs}.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J job_{filename}_{mapping}_{nprocs}
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

export OMPI_MCA_mpi_warn_on_fork=0
export OMPI_MCA_mpi_yield_when_idle=1

"""

filenames       = ['maxwell3d']
mappings        = [[['identity', True]]]
nnodes          = [1,2,2**2,2**3,2**4,2**5,2**6,2**7]
ntasks_per_node = 35
ncells          = [64,96,128,160,192,224]
degrees         = [2,3,4,5]

script_nc_d = 'srun python3 {filename}.py -n {nc} {nc} {nc} -d {d} {d} {d} {mapping}\n'

os.makedirs('mesh', exist_ok=True)
#os.system('python3 mesh/generate_meshes.py')

os.makedirs('results', exist_ok=True)

for mapping,f in zip(mappings,filenames):
    for ea in mapping:
        mapp = '-m ' + ea[0]+(' --a' if ea[1] else '')
        mapp_name = ea[0] + (('_'+'analytical') if ea[1] else '')

        for nn in nnodes:
            batch_script = batch_str.format(filename=f,mapping=mapp_name, nprocs=nn*ntasks_per_node,nnodes=nn,ntasks_per_node=ntasks_per_node, mem=150000,time_limit="6:00:00")
            for nc in ncells:
                for d in degrees:
                    batch_script += script_nc_d.format(filename=f, nc=nc, d=d, mapping=mapp)
                
                batch_script += '\n'

            filename = 'job_{filename}_{mapping}_{nprocs}.sh'.format(filename=f, mapping=mapp_name, nprocs=nn*ntasks_per_node)
            with open(filename,'w') as file_:
                file_.write(batch_script)

            os.system('sbatch {}'.format(filename))


    
    

        
    


