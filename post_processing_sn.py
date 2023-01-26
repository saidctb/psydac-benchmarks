import numpy as np


petiga_assembly_32 = np.array([[2.812678, 1.555927, 0.826999, 0.449456, 0.240059, 0.12946],
                     [94.812609, 48.436335, 24.678121, 12.678361, 6.68989, 3.44336],
                     [497.478098, 278.735268, 145.677756, 76.15136, 42.889353, 23.884055]])

petiga_assembly_40 = np.array([[5.517884, 2.987778, 1.573295, 0.846826, 0.450718, 0.233881],
                     [184.765686, 93.608816, 48.773257, 24.63235, 12.825916, 6.60611],
                     [963.656463, 541.673641, 286.931461, 145.815962, 84.127257, 46.483207]])

petiga_assembly = np.array([[[petiga_assembly_32, petiga_assembly_40]]])
petiga_assembly_scaling  = np.array([[[petiga_assembly[0,0,i1,i2,0]/petiga_assembly[0,0,i1,i2,:]/[1,2,4,8,16,32] for i2 in range(3)] for i1 in range(2)]])

petiga_dot_32 = np.array([[0.007299, 0.00418, 0.00216, 0.001216, 0.00068, 0.000286],
                 [0.109616, 0.056852, 0.029139, 0.015233, 0.008292, 0.00423],
                 [0.530624, 0.298215, 0.156179, 0.082217, 0.046765, 0.025973]])

petiga_dot_40 = np.array([[0.015556, 0.00834, 0.004335, 0.002327, 0.001415, 0.000661],
                 [0.21524, 0.109356, 0.057291, 0.029279, 0.015955, 0.008169],
                 [1.02939, 0.576643, 0.307164, 0.156558, 0.091616, 0.050182]])

petiga_dot = np.array([[[petiga_dot_32, petiga_dot_40]]])
petiga_dot_scaling  = np.array([[[petiga_dot[0,0,i1,i2,0]/petiga_dot[0,0,i1,i2,:]/[1,2,4,8,16,32] for i2 in range(3)] for i1 in range(2)]])
#######################################################################################################
import numpy as np

results_folder = 'results/'
problems = ['poisson_3d','time_harmonic_maxwell_3d','biharmonic_3d']
problems = ['biharmonic_3d']
mappings = [[['identity', True],['identity', False]],[['identity', True]],[['identity', False]]]
mappings = [[['identity', False]]]
ncells   = [32,40]
degrees  = [2,3,4]
number_of_nodes     = [1,1,1,1,1,1]
number_of_mpi_procs = np.array([1,2,4,8,16,32])
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
                        timmings_solve[i1,i2,i3,i4,i5]       = T['solve_time']
                        timmings_dot_p[i1,i2,i3,i4,i5]       = T['dot_product_time']
                    except:
                        timmings_bi_assembly[i1,i2,i3,i4,i5] = np.nan
                        timmings_dot_p[i1,i2,i3,i4,i5] = np.nan


                k = [j for j in range(len(number_of_mpi_procs)) if not np.isnan(timmings_bi_assembly[i1,i2,i3,i4,j])] + [0]
                nn = [np.nan]*k[0] + [2**(i-k[0]) for i in range(k[0],len(number_of_mpi_procs))]
                scaling_bi_assembly[i1,i2,i3,i4,:] = timmings_bi_assembly[i1,i2,i3,i4,k[0]]/timmings_bi_assembly[i1,i2,i3,i4,:]/nn
                scaling_dot_p[i1,i2,i3,i4,:]       = timmings_dot_p[i1,i2,i3,i4,k[0]]/timmings_dot_p[i1,i2,i3,i4,:]/nn

#########################################################################################################

number_of_nodes = [1,1,1,1,1,1]
number_of_mpi_procs = np.array([1,1,1,1,1,1])
number_of_threads = [1,2,4,8,16,32]

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

#########################################################################################################
from tabulate import tabulate
headers = [""] + [str(np*nt) for np,nt in zip(number_of_nodes, number_of_threads)]

for i1,p in enumerate(problems):
    for i2,mapping in enumerate(mappings[i1]):
        if all(np.isnan(v) for v in timmings_bi_assembly[i1,i2].flatten()):continue
        mapping = ('{} analytical mapping' if mapping[1] else '{} Nurbs mapping').format(mapping[0])
        print("="*45,"Timings of the Matrix Assembly of {} with the {}".format(p,mapping), "="*45)
        T = np.around(petiga_dot[i1,i2], decimals=5)
        newT = []
        for i3,nc in enumerate(ncells):
            for i4,d in enumerate(degrees):
                newT.append(["nc = {} ** 3 , p = {}".format(nc,d)] +  T[i3,i4].tolist())
            newT.append(["   "]*len(T[0]))
        
        print(tabulate(newT, headers=headers, tablefmt="grid"))
        print("\n")
raise
#from tabulate import tabulate
#headers = [""] + ['$p = {}$'.format(d) for d in degrees]
#paralle_ef = [[petiga_assembly_scaling]]
#titles = ['Matrix Assembly', 'Matrix Vector Product']

#for i1,p in enumerate(problems):
#    for i2,mapping in enumerate(mappings[i1]):
#        for paralle_ef_m,title in zip(paralle_ef, titles):
#            print("="*45,"Parallel Efficency of {}".format(title), "="*45)
#            T1 = np.around(paralle_ef[0][0], decimals=4)
#            newT1 = []
#            for i3,nc in enumerate(ncells):
#                newT1.append(["$ n_{{el}} = {}^3 $".format(nc)]+ ['{}%'.format(int(T1[i3,i4][-1]*10000)/100) for i4,d in enumerate(degrees)])

#            print(tabulate(newT1, headers=headers, tablefmt="latex"))
#            print("\n")
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
fnames  = ['matrix_assembly_biharmonic_strong_scaling_single_node', 'matrix_vector_product_biharmonic_strong_scaling_single_node','matrix_assembly_biharmonic_strong_scaling_single_node_multi_threading', 'matrix_vector_product_biharmonic_strong_scaling_single_node_multi_threading']
xaxist = [r'Number of processes or threads']*4
timings = [[timmings_bi_assembly[0,0], timmings_bi_assembly_mth[0,0], petiga_assembly], [timmings_dot_p[0,0], timmings_dot_p_mth[0,0], petiga_dot]]
number_of_mpi_procs = np.array([1,2,4,8,16,32])

for title,fname,timings_i,xlabel in zip(titles, fnames, timings,xaxist):
    fig = plt.figure(figsize=(10,15))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(number_of_mpi_procs,[5*np.nanmax(timings_i[1])/2**d for d in range(len(number_of_mpi_procs))], color='black', linestyle='dashed', label='perfect scaling')
    for nc in range(len(ncells)):
        for p in range(degrees[0],degrees[-1]+1):

            mask = np.isfinite(timings_i[0][nc,p-degrees[0]])
            line, = ax.plot(number_of_mpi_procs[mask], timings_i[0][nc,p-degrees[0]][mask], line_styles[nc],color=colors[p-degrees[0]])

            mask = np.isfinite(timings_i[1][nc,p-degrees[0]])
            line, = ax.plot(number_of_mpi_procs[mask], timings_i[1][nc,p-degrees[0]][mask], marker=markers[nc], linestyle='dashed',color=colors[p-degrees[0]])

            mask = np.isfinite(timings_i[2][nc][p-degrees[0]])
            line, = ax.plot(number_of_mpi_procs[mask], timings_i[2][nc][p-degrees[0]][mask], marker=markers[nc], linestyle='dotted', color=colors[p-degrees[0]])

        row = '$n_{{el}}={}^3$'.format(ncells[nc])
        line, = ax.plot(np.nan*number_of_mpi_procs[mask], np.nan*timings_i[0][nc,0][mask], line_styles[nc],color='k', label=row)

    for p in range(degrees[0],degrees[-1]+1):
        row = '$p={}$ (Psydac Pure MPI)'.format(p)
        line, = ax.plot(np.nan*number_of_mpi_procs[mask], np.nan*timings_i[0][0,p-degrees[0]][mask],color=colors[p-degrees[0]], label=row)
        row = '$p={}$ (Psydac OpenMP)'.format(p)
        line, = ax.plot(np.nan*number_of_mpi_procs[mask], np.nan*timings_i[0][0,p-degrees[0]][mask],linestyle='dashed', color=colors[p-degrees[0]], label=row)
        row = '$p={}$ (PetIGA Pure MPI)'.format(p)
        line, = ax.plot(np.nan*number_of_mpi_procs[mask], np.nan*timings_i[0][0,p-degrees[0]][mask],linestyle='dotted', color=colors[p-degrees[0]], label=row)

    box = ax.get_position()
#   ax.set_position([box.x0, box.y0, box.width * 0.3, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_xlabel( xlabel, rotation='horizontal' )
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

