import numpy as np

results_folder      = 'results/'
problems            = ['h1_projector_3d', 'hcurl_projector_3d', 'hdiv_projector_3d']
ncells              = [64,96,128,160,192,224]
degrees             = [2,3,4,5]
number_of_mpi_procs = [1*32,2*32,4*32,8*32,16*32,32*32,64*32,128*32]
number_of_threads   = 1

timmings_assembly = np.zeros((len(problems), len(ncells),len(degrees), len(number_of_mpi_procs)))
timmings_solver   = np.zeros((len(problems), len(ncells),len(degrees), len(number_of_mpi_procs)))
scaling_assembly  = timmings_assembly.copy()
scaling_solver    = timmings_solver.copy()

for i1,p in enumerate(problems):
    for i3,nc in enumerate(ncells):
        for i4,d in enumerate(degrees):
            for i5, mpi_p in enumerate(number_of_mpi_procs):
                names = (p,) + (nc,)*3 +(d,)*3 + (mpi_p,)
                filename = '_'.join([str(i) for i in names])
                try:
                    T = np.load(results_folder+filename+'.npy', allow_pickle=True)
                    T = T.item()
                    timmings_assembly[i1,i3,i4,i5] = T['matrix_assembly']
                    timmings_solver  [i1,i3,i4,i5] = T['solver']
                except:
                    timmings_assembly[i1,i3,i4,i5] = np.nan
                    timmings_solver  [i1,i3,i4,i5] = np.nan

            k = [j for j in range(len(number_of_mpi_procs)) if not np.isnan(timmings_assembly[i1,i3,i4,j])] + [0]
            scaling_assembly[i1,i3,i4,:] = timmings_assembly[i1,i3,i4,:]/timmings_assembly[i1,i3,i4,k[0]]
            scaling_solver  [i1,i3,i4,:] = timmings_solver  [i1,i3,i4,:]/timmings_solver  [i1,i3,i4,k[0]]

from tabulate import tabulate
headers = [""] + [str(np) for np in number_of_mpi_procs]

for i1,p in enumerate(problems):
    if all(np.isnan(v) for v in timmings_solver[i1].flatten()):continue
    print("="*45,"Timings of the kronecker solver of {}".format(p), "="*45)
    T = np.around(timmings_solver[i1], decimals=5)
    newT = []
    for i2,nc in enumerate(ncells):
        for i3,d in enumerate(degrees):
            newT.append(["nc = {} ** 3 , p = {}".format(nc,d)] +  T[i2,i3].tolist())
        newT.append(["   "]*len(T[0]))
    
    print(tabulate(newT, headers=headers, tablefmt="grid"))
    print("\n")
                
