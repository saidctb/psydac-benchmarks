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
        T = np.around(timmings_bi_assembly[i1,i2], decimals=5)
        newT = []
        for i2,nc in enumerate(ncells):
            for i3,d in enumerate(degrees):
                newT.append(["nc = {} ** 3 , p = {}".format(nc,d)] +  T[i2,i3].tolist())
            newT.append(["   "]*len(T[0]))
        
        print(tabulate(newT, headers=headers, tablefmt="grid"))
        print("\n")


#====================================================================================================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.legend_handler import HandlerLine2D

colors = np.linspace(0, 1, len(degrees))
colors = cm.rainbow(colors)
line_styles = ['>-','o-','s-','v-']

from itertools import product

titles = ['Matrix Assembly', 'Matrix Vector Product','Matrix Assembly', 'Matrix Vector Product']
fnames  = ['matrix_assembly_bih_sn', 'matrix_vector_product_bih_sn','matrix_assembly_mth_bih_sn', 'matrix_vector_product_mth_bih_sn']
xaxist = [r'number of mpi procs', r'number of mpi procs',r'number of threads', r'number of threads']
timings = [timmings_bi_assembly[0,0], timmings_dot_p[0,0],timmings_bi_assembly_mth[0,0], timmings_dot_p_mth[0,0]]
number_of_mpi_procs = np.array([1,2,4,8,16,32])

for title,fname,timings_i,xlabel in zip(titles, fnames, timings,xaxist):
    fig = plt.figure(figsize=(10,15))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(number_of_mpi_procs,[5*np.nanmax(timings_i)/2**d for d in range(len(number_of_mpi_procs))], color='black', linestyle='dashed', label='perfect scaling')
    for nc in range(len(ncells)):
        for p in range(degrees[0],degrees[-1]+1):

            mask = np.isfinite(timings_i[nc,p-degrees[0]])
            line, = ax.plot(number_of_mpi_procs[mask], timings_i[nc,p-degrees[0]][mask], line_styles[nc],color=colors[p-degrees[0]])

        row = '$n_{{el}}={}^3$'.format(ncells[nc])
        line, = ax.plot(np.nan*number_of_mpi_procs[mask], np.nan*timings_i[nc,degrees[0]][mask], line_styles[nc],color='k', label=row)

    for p in range(degrees[0],degrees[-1]+1):
        row = '$p={}$'.format(p)
        line, = ax.plot(np.nan*number_of_mpi_procs[mask], np.nan*timings_i[0,degrees[0]][mask],color=colors[p-degrees[0]], label=row)


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

    fig.savefig(fname)
