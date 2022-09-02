import numpy as np

results_folder = 'results/'
problems = ['poisson_3d','vector_poisson_3d','time_harmonic_maxwell_3d']
mappings = [[['identity', True],['identity', False], ['quarter_annulus', True], ['quarter_annulus', False]],[['identity', True]],[['identity', True]]]
ncells   = [10,96,128,160,192,256]
degrees  = [2,3,4,5]
number_of_mpi_procs = [1,1*32,2*32,4*32,8*32,16*32,32*32,64*32,128*32]
number_of_threads = 1

timmings_bi_assembly   = np.zeros((len(problems), max(len(mapping) for mapping in mappings), len(ncells),len(degrees), len(number_of_mpi_procs)))
timmings_dot_p         = np.zeros((len(problems), max(len(mapping) for mapping in mappings), len(ncells),len(degrees), len(number_of_mpi_procs)))
scaling_bi_assembly    = timmings_bi_assembly.copy()
scaling_dot_p          = timmings_dot_p.copy()

for i1,p in enumerate(problems):
    for i2,mapping in enumerate(mappings[i1]):
        for i3,nc in enumerate(ncells):
            for i4,d in enumerate(degrees):
                for i5, mpi_p in enumerate(number_of_mpi_procs):
                    names = (p, mapping[0])+ (('geof',) if not mapping[1] else ()) + (nc,)*3 +(d,)*3 + (mpi_p, number_of_threads)
                    filename = '_'.join([str(i) for i in names])
                    try:
                        T = np.load(results_folder+filename+'.npy', allow_pickle=True)
                        T = T.item()
                        timmings_bi_assembly[i1,i2,i3,i4,i5] = T['bilinear_form_assembly_time']
                        timmings_dot_p[i1,i2,i3,i4,i5]       = T['dot_product_time']
                    except:
                        timmings_bi_assembly[i1,i2,i3,i4,i5] = np.nan
                        timmings_dot_p[i1,i2,i3,i4,i5] = np.nan

                k = [j for j in range(len(number_of_mpi_procs)) if not np.isnan(timmings_bi_assembly[i1,i2,i3,i4,j])] + [0]
                scaling_bi_assembly[i1,i2,i3,i4,:] = timmings_bi_assembly[i1,i2,i3,i4,:]/timmings_bi_assembly[i1,i2,i3,i4,k[0]]
                scaling_dot_p[i1,i2,i3,i4,:]       = timmings_dot_p[i1,i2,i3,i4,:]/timmings_dot_p[i1,i2,i3,i4,k[0]]

from tabulate import tabulate
headers = [""] + [str(np) for np in number_of_mpi_procs]

for i1,p in enumerate(problems):
    for i2,mapping in enumerate(mappings[i1]):
        if all(np.isnan(v) for v in timmings_bi_assembly[i1,i2].flatten()):continue
        mapping = ('{} analytical mapping' if mapping[1] else '{} Nurbs mapping').format(mapping[0])
        print("="*45,"Timings of the Matrix Assembly of {} with the {}".format(p,mapping), "="*45)
        T = np.around(timmings_bi_assembly[i1,i2], decimals=5)
        newT = []
        for i2,nc in enumerate(ncells):
            for i3,d in enumerate(degrees):
                newT.append(["nc = {} ** 3 , p = {}".format(nc,d)] +  T[i2,i3].tolist())
            newT.append(["   "]*len(T[0]))
        
        print(tabulate(newT, headers=headers, tablefmt="grid"))
        print("\n")
                
