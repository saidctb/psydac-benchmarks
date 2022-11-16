import numpy as np

results_folder = 'results/'
ncells   = [50,75,100,125,150, 175,200]
degrees  = [2,3]
number_of_mpi_procs = [1*35,2*35,4*35,8*35,16*35,32*35,64*35,128*35]
number_of_threads   = 1
filename            = 'maxwell3d'

timmings_bi_assembly = np.zeros((len(ncells),len(degrees), len(number_of_mpi_procs)))
timmings_dot_p       = np.zeros((len(ncells),len(degrees), len(number_of_mpi_procs)))
scaling_bi_assembly  = timmings_bi_assembly.copy()
scaling_dot_p        = timmings_dot_p.copy()

for i1,nc in enumerate(ncells):
    for i2,d in enumerate(degrees):
        for i3, mpi_p in enumerate(number_of_mpi_procs):
            names = (p,) + (nc,)*3 +(d,)*3 + (mpi_p, number_of_threads)
            filename = '_'.join([str(i) for i in names])
            try:
                T = np.load(results_folder+filename+'.npy', allow_pickle=True)
                T = T.item()
                timmings_bi_assembly[i1,i2,i3] = T['bilinear_form_assembly_time']
                timmings_dot_p[i1,i2,i3]       = T['dot_product_time']
            except:
                timmings_bi_assembly[i1,i2,i3] = np.nan
                timmings_dot_p[i1,i2,i3] = np.nan

        k = [j for j in range(len(number_of_mpi_procs)) if not np.isnan(timmings_bi_assembly[i1,i2,j])] + [0]
        scaling_bi_assembly[i1,i2,:] = timmings_bi_assembly[i1,i2,:]/timmings_bi_assembly[i1,i2,k[0]]
        scaling_dot_p[i1,i2,:]       = timmings_dot_p[i1,i2,:]/timmings_dot_p[i1,i2,k[0]]

#-----------------------------------------------------------------------------------------------------------------------
from tabulate import tabulate
headers = [""] + [str(np) for np in number_of_mpi_procs]

if not all(np.isnan(v) for v in timmings_bi_assembly.flatten()):
    print("="*45,"Timings of the Matrix Assembly of Time dependent Maxwell", "="*45)
    T = np.around(timmings_bi_assembly, decimals=5)
    newT = []
    for i1,nc in enumerate(ncells):
        for i2,d in enumerate(degrees):
            newT.append(["nc = {} ** 3 , p = {}".format(nc,d)] +  T[i1,i2].tolist())
        newT.append(["   "]*len(T[0]))

    print(tabulate(newT, headers=headers, tablefmt="grid"))
    print("\n")
    
if not all(np.isnan(v) for v in timmings_dot_p.flatten()):
    print("="*45,"Timings of the Matrix vector dot product of the Hcurl Mass Matrix", "="*45)
    T = np.around(timmings_dot_p, decimals=5)
    newT = []
    for i1,nc in enumerate(ncells):
        for i2,d in enumerate(degrees):
            newT.append(["nc = {} ** 3 , p = {}".format(nc,d)] +  T[i1,i2].tolist())
        newT.append(["   "]*len(T[0]))

    print(tabulate(newT, headers=headers, tablefmt="grid"))
    print("\n")
            
