import numpy as np
import os

problems = ['poisson_3d','vector_poisson_3d','time_harmonic_maxwell_3d']
mappings = [[['identity', True]],[['identity', True]],[['identity', True]]]
ncells   = [32,64,96,128,160,192,256]
degrees  = [2,3,4,5]
number_of_mpi_procs = [1*32,2*32,4*32,8*32,16*32,32*32,64*32,128*32]
number_of_threads = 1

for i1,p in enumerate(problems):
    for i2,mapping in enumerate(mappings[i1]):
        for i3,nc in enumerate(ncells):
            for i4,d in enumerate(degrees):
                for i5, mpi_p in enumerate(number_of_mpi_procs):
                    names = (p,) + (nc,)*3 +(d,)*3 + (mpi_p, number_of_threads)
                    oldfilename = '_'.join([str(i) for i in names])
                    names = (p, mapping[0])+ (('geof',) if not mapping[1] else ()) + (nc,)*3 +(d,)*3 + (mpi_p, number_of_threads)
                    newfilename = '_'.join([str(i) for i in names])
                    
                    try:
                        T = np.load(oldfilename+'.npy', allow_pickle=True)
                        T = T.item()
                        np.save(newfilename, T)
                        os.remove(oldfilename+'.npy')
                    except:
                        pass
