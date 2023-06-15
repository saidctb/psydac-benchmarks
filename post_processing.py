import numpy as np
import re
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})

results_folder = 'results/'
petiga_results_folder = 'petiga2/'
get_digits = lambda s:tuple(float(i) for i in re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", s)[-2:])

#problems = ['poisson_3d','time_harmonic_maxwell_3d','biharmonic_3d']
problems = ['biharmonic_3d']
#mappings = [[['identity', True],['identity', False]],[['identity', True]],[['identity', False]]]
mappings = [[['quarter_annulus', False]]]
ncells   = [128,160]
degrees  = [2,3,4]
number_of_nodes     = np.array([1,2,4,8,16,32,64])
number_of_mpi_procs = number_of_nodes*32
number_of_threads   = 1

timmings_bi_assembly   = np.zeros((len(problems), max(len(mapping) for mapping in mappings), len(ncells),len(degrees), len(number_of_mpi_procs)))
timmings_dot_p         = np.zeros((len(problems), max(len(mapping) for mapping in mappings), len(ncells),len(degrees), len(number_of_mpi_procs)))
timmings_solve         = np.zeros((len(problems), max(len(mapping) for mapping in mappings), len(ncells),len(degrees), len(number_of_mpi_procs)))
scaling_bi_assembly    = timmings_bi_assembly.copy()
scaling_dot_p          = timmings_dot_p.copy()
scaling_solve          = timmings_solve.copy()

for i1,p in enumerate(problems):
    for i2,mapping in enumerate(mappings[i1]):
        for i3,nc in enumerate(ncells):
            for i4,d in enumerate(degrees):
                for i5, mpi_p in enumerate(number_of_mpi_procs):
                    names = (p,)+ (('geof',) if not mapping[1] else ()) + (nc,)*3 +(d,)*3 + (mpi_p, number_of_threads)
                    filename = '_'.join([str(i) for i in names])
                    try:
                        T = np.load(results_folder+filename+'.npy', allow_pickle=True)
                        T = T.item()
                        timmings_bi_assembly[i1,i2,i3,i4,i5] = min(T['bilinear_form_assembly_time'],T.get('blinear_form_assembly_time2', T['bilinear_form_assembly_time']))
                        timmings_dot_p[i1,i2,i3,i4,i5]       = T['dot_product_time']
                    except:
                        timmings_bi_assembly[i1,i2,i3,i4,i5] = np.nan
                        timmings_dot_p[i1,i2,i3,i4,i5] = np.nan


                k = [j for j in range(len(number_of_mpi_procs)) if not np.isnan(timmings_bi_assembly[i1,i2,i3,i4,j])] + [0]
                nn = [np.nan]*k[0] + [2**(i-k[0]) for i in range(k[0],len(number_of_mpi_procs))]
                scaling_bi_assembly[i1,i2,i3,i4,:] = timmings_bi_assembly[i1,i2,i3,i4,k[0]]/timmings_bi_assembly[i1,i2,i3,i4,:]/nn
                scaling_dot_p[i1,i2,i3,i4,:]       = timmings_dot_p[i1,i2,i3,i4,k[0]]/timmings_dot_p[i1,i2,i3,i4,:]/nn

#########################################################################################################

timmings_bi_assembly_pet = np.zeros((len(ncells),len(degrees), len(number_of_mpi_procs)))
timmings_dot_p_pet       = np.zeros((len(ncells),len(degrees), len(number_of_mpi_procs)))
scaling_bi_assembly_pet  = timmings_bi_assembly_pet.copy()
scaling_dot_p_pet        = timmings_dot_p_pet.copy()
for i1,nc in enumerate(ncells):
    for i2,d in enumerate(degrees):
        for i3, mpi_p in enumerate(number_of_mpi_procs):
            names = ('bi_quarter_annulus_3d', nc, d, mpi_p)
            filename = '_'.join([str(i) for i in names])
            try:
                fname = petiga_results_folder+filename+'.out'
                f = open(fname, 'r')
                digits = get_digits(f.read())
                timmings_bi_assembly_pet[i1,i2,i3] = digits[0]
                timmings_dot_p_pet[i1,i2,i3]       = digits[1]
            except:
                raise
                timmings_bi_assembly_pet[i1,i2,i3] = np.nan
                timmings_dot_p_pet[i1,i2,i3,i4,i5] = np.nan


        k = [j for j in range(len(number_of_mpi_procs)) if not np.isnan(timmings_bi_assembly_pet[i1,i2,j])] + [0]
        nn = [np.nan]*k[0] + [2**(i-k[0]) for i in range(k[0],len(number_of_mpi_procs))]
        scaling_bi_assembly_pet[i1,i2,:] = timmings_bi_assembly_pet[i1,i2,k[0]]/timmings_bi_assembly_pet[i1,i2,:]/nn
        scaling_dot_p_pet[i1,i2,:]       = timmings_dot_p_pet[i1,i2,k[0]]/timmings_dot_p_pet[i1,i2,:]/nn

number_of_mpi_procs = number_of_nodes*1
number_of_threads = [32]*len(number_of_nodes)

timmings_bi_assembly_mth  = np.zeros((len(problems), max(len(mapping) for mapping in mappings), len(ncells),len(degrees), len(number_of_mpi_procs)))
timmings_dot_p_mth        = np.zeros((len(problems), max(len(mapping) for mapping in mappings), len(ncells),len(degrees), len(number_of_mpi_procs)))
timmings_solve_mth        = np.zeros((len(problems), max(len(mapping) for mapping in mappings), len(ncells),len(degrees), len(number_of_mpi_procs)))
scaling_bi_assembly_mth   = timmings_bi_assembly_mth.copy()
scaling_dot_p_mth         = timmings_dot_p_mth.copy()
scaling_solve_mth         = timmings_solve.copy()

for i1,p in enumerate(problems):
    for i2,mapping in enumerate(mappings[i1]):
        for i3,nc in enumerate(ncells):
            for i4,d in enumerate(degrees):
                for i5, mpi_p in enumerate(number_of_mpi_procs):
                    names = (p,)+ (('geof',) if not mapping[1] else ()) + (nc,)*3 +(d,)*3 + (mpi_p, number_of_threads[i5])
                    filename = '_'.join([str(i) for i in names])
                    try:
                        T = np.load(results_folder+filename+'.npy', allow_pickle=True)
                        T = T.item()
                        timmings_bi_assembly_mth[i1,i2,i3,i4,i5] = min(T['bilinear_form_assembly_time'],T.get('blinear_form_assembly_time2', T['bilinear_form_assembly_time']))
                        timmings_dot_p_mth[i1,i2,i3,i4,i5]       = T['dot_product_time']
                    except:
                        timmings_bi_assembly_mth[i1,i2,i3,i4,i5] = np.nan
                        timmings_dot_p_mth[i1,i2,i3,i4,i5] = np.nan


                k = [j for j in range(len(number_of_mpi_procs)) if not np.isnan(timmings_bi_assembly_mth[i1,i2,i3,i4,j])] + [0]
                nn = [np.nan]*k[0] + [2**(i-k[0]) for i in range(k[0],len(number_of_mpi_procs))]
                scaling_bi_assembly_mth[i1,i2,i3,i4,:] = timmings_bi_assembly_mth[i1,i2,i3,i4,k[0]]/timmings_bi_assembly_mth[i1,i2,i3,i4,:]/nn
                scaling_dot_p_mth[i1,i2,i3,i4,:]       = timmings_dot_p_mth[i1,i2,i3,i4,k[0]]/timmings_dot_p_mth[i1,i2,i3,i4,:]/nn

##########################################################################################################
#from tabulate import tabulate
#headers = [""] + ['$p = {}$'.format(d) for d in degrees]
#paralle_ef = [[scaling_bi_assembly,scaling_bi_assembly_mth]]
#titles = ['Matrix Assembly', 'Matrix Vector Product']

#for i1,p in enumerate(problems):
#    for i2,mapping in enumerate(mappings[i1]):
#        for paralle_ef_m,title in zip(paralle_ef, titles):
#            print("="*45,"Parallel Efficency of {}".format(title), "="*45)
#            T1 = np.around(paralle_ef_m[0][i1,i2], decimals=4)
#            T2 = np.around(paralle_ef_m[1][i1,i2], decimals=4)
#            newT1 = []
#            newT2 = []
#            for i3,nc in enumerate(ncells):
#                newT1.append(["$ n_{{el}} = {}^3 $".format(nc)]+ ['{}%'.format(int(T1[i3,i4][-1]*10000)/100) for i4,d in enumerate(degrees)])
#                newT2.append(["$ n_{{el}} = {}^3 $".format(nc)]+ ['{}%'.format(int(T2[i3,i4][-1]*10000)/100) for i4,d in enumerate(degrees)])

#            print(tabulate(newT1, headers=headers, tablefmt="latex"))
#            print(tabulate(newT2, headers=headers, tablefmt="latex"))
#            print("\n")
#raise
#from tabulate import tabulate
#headers = [""] + [str(nn) for nn in number_of_nodes]

#timmings_dot_p_mth[0,0] = np.array([[[0.00943, 0.00495, 0.00266, 0.00164, 0.00091, 0.00054, 0.00042],
#        [0.01851, 0.00966, 0.00513, 0.00294, 0.00163, 0.0009 , 0.00059],
#        [0.04284, 0.02241, 0.01189, 0.00643, 0.00349, 0.00217, 0.00111]],

#       [[0.01827, 0.0094 , 0.00497, 0.00292, 0.00159, 0.00093, 0.00075],
#        [0.03596, 0.01847, 0.00964, 0.00532, 0.00288, 0.00161, 0.001  ],
#        [0.08815, 0.04648, 0.02309, 0.01191, 0.00634, 0.00332, 0.00187]]])
#for i1,p in enumerate(problems):
#    for i2,mapping in enumerate(mappings[i1]):
#        if all(np.isnan(v) for v in timmings_bi_assembly[i1,i2].flatten()):continue
#        mapping = ('{} analytical mapping' if mapping[1] else '{} Nurbs mapping').format(mapping[0])
#        print("="*45,"Timings of the Matrix Assembly of {} with the {}".format(p,mapping), "="*45)
#        T = np.around(timmings_dot_p_mth[i1,i2], decimals=5)
#        newT = []
#        for i3,nc in enumerate(ncells):
#            for i4,d in enumerate(degrees):
#                newT.append(["nc = {} ** 3 , p = {}".format(nc,d)] +  T[i3,i4].tolist())
#            newT.append(["   "]*len(T[0]))
#        
#        print(tabulate(newT, headers=headers, tablefmt="grid"))
#        print("\n")
#raise
#====================================================================================================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.legend_handler import HandlerLine2D

colors = np.linspace(0, 1, len(degrees))
colors = cm.rainbow(colors)
line_styles = ['>-','o-','s-','v-']
markers = ['>','o','s','v']

from itertools import product

titles = ['Matrix Assembly', 'Matrix Vector Product','Matrix Assembly', 'Matrix Vector Product']
fnames = ['matrix_assembly_biharmonic_strong_scaling', 'matrix_vector_product_biharmonic_strong_scaling','matrix_assembly_biharmonic_strong_scaling_multi_threading', 'matrix_vector_product_biharmonic_strong_scaling_multi_threading']
xaxist = [r'number of nodes', r'number of nodes',r'number of nodes',r'number of nodes']
timings = [[timmings_bi_assembly[0,0], timmings_bi_assembly_mth[0,0]], [timmings_dot_p[0,0], timmings_dot_p_mth[0,0]]]


#titles = ['Matrix Vector Product']
#fnames = ['matrix_vector_product_mth_bih_ss']
#xaxist = [r'number of nodes']
#timings = [timmings_dot_p_mth[0,0]]
number_of_nodes = np.array([1,2,4,8,16,32,64])

for title,fname,timings_i, xlabel in zip(titles, fnames, timings,xaxist):
    fig = plt.figure(figsize=(15,8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(number_of_nodes,[5*np.nanmax(timings_i[1])/2**d for d in range(len(number_of_mpi_procs))], color='black', linestyle='dashed', label='perfect scaling')
    for nc in range(len(ncells)):
        for p in range(degrees[0],degrees[-1]+1):

            mask = np.isfinite(timings_i[0][nc,p-degrees[0]])
            line, = ax.plot(number_of_nodes[mask], timings_i[0][nc,p-degrees[0]][mask], line_styles[nc],color=colors[p-degrees[0]])

            mask = np.isfinite(timings_i[1][nc,p-degrees[0]])
            line, = ax.plot(number_of_nodes[mask], timings_i[1][nc,p-degrees[0]][mask], marker=markers[nc], linestyle='dashed', color=colors[p-degrees[0]])

        row = '$n_{{el}}={}^3$'.format(ncells[nc])
        line, = ax.plot(np.nan*number_of_nodes[mask], np.nan*timings_i[0][nc,0][mask], line_styles[nc],color='k', label=row)

    for p in range(degrees[0],degrees[-1]+1):
        row = '$p={}$ (Pure MPI)'.format(p)
        line, = ax.plot(np.nan*number_of_nodes[mask], np.nan*timings_i[0][0,p-degrees[0]][mask],color=colors[p-degrees[0]], label=row)
        row = '$p={}$ (MPI+OpenMP)'.format(p)
        line, = ax.plot(np.nan*number_of_nodes[mask], np.nan*timings_i[0][0,p-degrees[0]][mask],linestyle='dashed', color=colors[p-degrees[0]], label=row)


    box = ax.get_position()
#   ax.set_position([box.x0, box.y0, box.width * 0.3, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_xlabel( r'number of nodes', rotation='horizontal' )
    ax.set_ylabel( r'time [s]' )
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xticks([])
    ax.set_xticks(number_of_mpi_procs)
    ax.set_xticklabels([str(d) for d in number_of_mpi_procs])
    ax.grid(True)
#    ax.title.set_text(title)
    fig.tight_layout(rect=[0, 0.05, 1, 1])

    fig.savefig("images/"+fname)

