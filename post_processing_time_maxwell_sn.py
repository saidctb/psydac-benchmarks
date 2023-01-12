import numpy as np

results_folder = 'results/'
filename       = 'maxwell3d'

number_of_mpi_procs = 7
ncells          = [52,64]
degrees         = [3,4,5]
number_of_threads = np.array([1,2,4,8,16])
timmings_bi_assembly     = np.zeros((len(ncells),len(degrees), len(number_of_threads)))
timmings_time_integrator = np.zeros((len(ncells),len(degrees), len(number_of_threads)))
timmings_dot_p           = np.zeros((len(ncells),len(degrees), len(number_of_threads)))
scaling_bi_assembly      = timmings_bi_assembly.copy()
scaling_time_integrator  = timmings_time_integrator.copy()
scaling_dot_p            = timmings_dot_p.copy()
for i1,nc in enumerate(ncells):
    for i2,d in enumerate(degrees):
        for i3, thr_p in enumerate(number_of_threads):
            names = ('maxwell_3d',) + (nc,)*3 +(d,)*3 + (number_of_mpi_procs, thr_p)
            filename = '_'.join([str(i) for i in names])
            try:
                T = np.load(results_folder+filename+'.npy', allow_pickle=True)
                T = T.item()
                timmings_bi_assembly[i1,i2,i3]     = T['bilinear_form_assembly_time']
                timmings_time_integrator[i1,i2,i3] = T['solution']
                timmings_dot_p[i1,i2,i3]           = T['dot_product_time']
            except:
                timmings_bi_assembly[i1,i2,i3] = np.nan
                timmings_time_integrator[i1,i2,i3] = np.nan
                timmings_dot_p[i1,i2,i3] = np.nan


        k = [j for j in range(len(number_of_threads)) if not np.isnan(timmings_bi_assembly[i1,i2,j])] + [0]
        nn = [np.nan]*k[0] + [2**(i-k[0]) for i in range(k[0],len(number_of_threads))]
        scaling_bi_assembly[i1,i2,:]     = timmings_bi_assembly[i1,i2,k[0]]/timmings_bi_assembly[i1,i2,:]/nn
        scaling_time_integrator[i1,i2,:] = timmings_time_integrator[i1,i2,k[0]]/timmings_time_integrator[i1,i2,:]/nn
        scaling_dot_p[i1,i2,:]           = timmings_dot_p[i1,i2,k[0]]/timmings_dot_p[i1,i2,:]/nn
#====================================================================================================
#from tabulate import tabulate
#headers = [""] + [str(nt) for nt in number_of_threads]

#if not all(np.isnan(v) for v in scaling_dot_p.flatten()):
#    print("="*45,"Timings of the Matrix Assembly of Time dependent Maxwell", "="*45)
#    T = np.around(scaling_bi_assembly, decimals=2)
#    newT = []
#    for i1,nc in enumerate(ncells):
#        for i2,d in enumerate(degrees):
#            newT.append(["nc = {} ** 3 , p = {}".format(nc,d)] +  T[i1,i2].tolist())
#        newT.append(["   "]*len(T[0]))

#    print(tabulate(newT, headers=headers, tablefmt="grid"))
#    print("\n")
#raise
##====================================================================================================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.legend_handler import HandlerLine2D

from itertools import product

colors = np.linspace(0, 1, len(degrees))
#np.random.shuffle(colors)
colors = cm.rainbow(colors)

titles = ['Matrix Assembly', 'Matrix Vector Product']
fnames = ['matrix_assembly_time_maxwell_single_node_multi_threading', 'matrix_vector_product_time_maxwell_single_node_multi_threading']
xaxist = [r'number of threads', r'number of threads']
timings = [timmings_bi_assembly, timmings_dot_p]

line_styles = ['>-','o-','s-','v-']
for title,fname,timings_i,xlabel in zip(titles, fnames, timings,xaxist):
    fig = plt.figure(figsize=(10,17))
#    fig.suptitle(title)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(number_of_threads,[5*np.nanmax(timings_i)/2**d for d in range(len(number_of_threads))], color='black', linestyle='dashed', label='perfect scaling')
    for nc in range(len(ncells)):
        for p in range(degrees[0],degrees[-1]+1):
            mask = np.isfinite(timings_i[nc,p-degrees[0]])
            line, = ax.plot(number_of_threads[mask], timings_i[nc,p-degrees[0]][mask], line_styles[nc],color=colors[p-degrees[0]])

        row = '$n_{{el}}={}^3$'.format(ncells[nc])
        line, = ax.plot(np.nan*number_of_threads[mask], np.nan*timings_i[nc,p-degrees[0]][mask], line_styles[nc],color='k', label=row)
            #ax.plot(number_of_threads,timings[nc,p-degrees[0]],'s', color=colors[p-degrees[0]])

    for p in range(degrees[0],degrees[-1]+1):
        row = '$p={}$'.format(p)
        line, = ax.plot(np.nan*number_of_threads[mask], np.nan*timings_i[0,p-degrees[0]][mask],color=colors[p-degrees[0]], label=row)

        # Shrink current axis by 20%
        box = ax.get_position()
#        ax.set_position([box.x0, box.y0, box.width * 0.3, box.height])

        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax.set_xlabel( xaxist, rotation='horizontal' )
        ax.set_ylabel( r'time [s]' )
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xticks([])
        ax.set_xticks(number_of_threads)
        ax.set_xticklabels([str(d) for d in number_of_threads])
        ax.grid(True)
#        ax.title.set_text('$ncells={}^3$'.format(ncells[nc]))
    fig.tight_layout()
    fig.savefig("images/"+fname)
