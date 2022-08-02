import numpy as np

results_folder = 'results/'
problems = ['poisson_3d','vector_poisson_3d','time_harmonic_maxwell_3d']
ncells   = [32,64,128,256]
degrees  = [2,3,4,5]
number_of_mpi_procs = [1*32,2*32,4*32,8*32,16*32,32*32,64*32,128*32]
number_of_threads = 1

timmings_bi_assembly   = np.zeros((len(problems), len(ncells),len(degrees), len(number_of_mpi_procs)))
timmings_dot_p         = np.zeros((len(problems), len(ncells),len(degrees), len(number_of_mpi_procs)))
scaling_bi_assembly    = timmings_bi_assembly.copy()
scaling_dot_p          = timmings_dot_p.copy()

for i1,p in enumerate(problems):
    for i2,nc in enumerate(ncells):
        for i3,d in enumerate(degrees):
            for i4, mpi_p in enumerate(number_of_mpi_procs):
                names = (p,) + (nc,)*3 +(d,)*3 + (mpi_p, number_of_threads)
                filename = '_'.join([str(i) for i in names])
                try:
                    T = np.load(results_folder+filename+'.npy', allow_pickle=True)
                    T = T.item()
                    timmings_bi_assembly[i1,i2,i3,i4] = T['bilinear_form_assembly_time']
                    timmings_dot_p[i1,i2,i3,i4]       = T['dot_product_time']
                except:
                    timmings_bi_assembly[i1,i2,i3,i4] = np.nan
                    timmings_dot_p[i1,i2,i3,i4] = np.nan

            k = [j for j in range(len(number_of_mpi_procs)) if not np.isnan(timmings_bi_assembly[i1,i2,i3,j])] + [0]
            scaling_bi_assembly[i1,i2,i3,:] = timmings_bi_assembly[i1,i2,i3,:]/timmings_bi_assembly[i1,i2,i3,k[0]]
            scaling_dot_p[i1,i2,i3,:]       = timmings_dot_p[i1,i2,i3,:]/timmings_dot_p[i1,i2,i3,k[0]]

from tabulate import tabulate
headers = [""] + [str(np) for np in number_of_mpi_procs]

for i1,p in enumerate(problems):
    print("="*45,"Timings of the Matrix Assembly of {}".format(p), "="*45)
    T = np.around(timmings_bi_assembly[i1], decimals=5).reshape((len(ncells)*len(degrees), len(number_of_mpi_procs))).tolist()
    newT = []
    for i2,nc in enumerate(ncells):
        for i3,d in enumerate(degrees):
            newT.append(["nc ={} ** 3 , p = {} ** 3".format(nc,d)] +  T[i2*len(ncells)+i3])
        newT.append(["   "]*len(T[0]))
    
    print(tabulate(newT, headers=headers, tablefmt="grid"))
    print("\n")
                
