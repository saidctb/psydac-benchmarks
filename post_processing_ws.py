import numpy as np

results_folder = 'results/'
f       = 'biharmonic_3d'

nnodes           = [1,1,1,1,1,1,2,4,8,16,32,64]
ntasks_per_node  = [1,2,4,8,16,32,32,32,32,32,32,32]
nthreads         = [1,1,1,1,1,1,1,1,1,1,1,1]

ncells0 = [[8,8,8],[16,8,8],[16,16,8],[16,16,16],[32,16,16],[32,32,16],[32,32,32],[64,32,32],[64,64,32],[64,64,64],[128,64,64],[128,128,64]]
ncells1 = [[10,10,10],[20,10,10],[20,20,10],[20,20,20],[40,20,20],[40,40,20],[40,40,40],[80,40,40],[80,80,40],[80,80,80],[160,80,80],[160,160,80]]
ncells = [ncells0,ncells1]
degrees = [2,3,4]

timmings_bi_assembly     = np.zeros((len(ncells),len(degrees), len(nnodes)))
timmings_dot_p           = np.zeros((len(ncells),len(degrees), len(nnodes)))

for i1,nc in enumerate(ncells):
    for i2,d in enumerate(degrees):
        for i3, (nn,nt,nth) in enumerate(zip(nnodes, ntasks_per_node, nthreads)):
            mpi_p = nn*nt
            names = (f,'geof') + (nc[i3][0],nc[i3][1],nc[i3][2]) +(d,)*3 + (mpi_p, nth)
            filename = '_'.join([str(i) for i in names])
            try:
                T = np.load(results_folder+filename+'.npy', allow_pickle=True)
                T = T.item()
                timmings_bi_assembly[i1,i2,i3]     = min(T['bilinear_form_assembly_time'],T.get('blinear_form_assembly_time2', T['bilinear_form_assembly_time']))
                timmings_dot_p[i1,i2,i3]           = T['dot_product_time']
            except:
                timmings_bi_assembly[i1,i2,i3] = np.nan
                timmings_solve_time[i1,i2,i3] = np.nan
                timmings_dot_p[i1,i2,i3] = np.nan

ntasks_per_node = [1,1,1,1,1,2,2,2,2,2,2,2]
nthreads        = [1,2,4,8,16,16,16,16,16,16,16,16]

timmings_bi_assembly_mth     = np.zeros((len(ncells),len(degrees), len(nnodes)))
timmings_dot_p_mth           = np.zeros((len(ncells),len(degrees), len(nnodes)))

for i1,nc in enumerate(ncells):
    for i2,d in enumerate(degrees):
        for i3, (nn,nt,nth) in enumerate(zip(nnodes, ntasks_per_node, nthreads)):
            mpi_p = nn*nt
            names = (f,'geof') + (nc[i3][0],nc[i3][1],nc[i3][2]) +(d,)*3 + (mpi_p, nth)
            filename = '_'.join([str(i) for i in names])

            try:
                T = np.load(results_folder+filename+'.npy', allow_pickle=True)
                T = T.item()
                timmings_bi_assembly_mth[i1,i2,i3]     = min(T['bilinear_form_assembly_time'],T.get('blinear_form_assembly_time2', T['bilinear_form_assembly_time']))
                timmings_dot_p_mth[i1,i2,i3]           = T['dot_product_time']
            except:
                timmings_bi_assembly_mth[i1,i2,i3] = np.nan
                timmings_dot_p_mth[i1,i2,i3] = np.nan

##-----------------------------------------------------------------------------------------------------------------------
#from tabulate import tabulate
#headers = [""] + [str(nn*nt*nth) for nn,nt,nth in zip(nnodes, ntasks_per_node, nthreads)]

#if not all(np.isnan(v) for v in timmings_bi_assembly.flatten()):
#    print("="*45,"Timings of the Matrix Assembly", "="*45)
#    T = np.around(timmings_dot_p_mth, decimals=5)
#    newT = []
#    for i1,nc in enumerate(ncells):
#        for i2,d in enumerate(degrees):
#            newT.append(["nc = {} ** 3 , p = {}".format(nc[0][0],d)] +  T[i1,i2].tolist())
#        newT.append(["   "]*len(T[0]))

#    print(tabulate(newT, headers=headers, tablefmt="grid"))
#    print("\n")
#raise
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
fnames = ['matrix_assembly_biharmonic_weak_scaling', 'matrix_vector_product_biharmonic_weak_scaling','matrix_assembly_biharmonic_weak_scaling_multi_threading', 'matrix_vector_product_biharmonic_weak_scaling_multi_threading']
xaxist = [r'number of mpi procs', r'number of mpi procs',r'number of threads',r'number of threads']
timings = [timmings_bi_assembly, timmings_dot_p, timmings_bi_assembly_mth, timmings_dot_p_mth]
nthreads = np.array([nn*nt*nth for nn,nt,nth in zip(nnodes, ntasks_per_node, nthreads)])
for title,fname,timings_i,xlabel in zip(titles, fnames, timings,xaxist):
    fig = plt.figure(figsize=(10,15))
    ax = fig.add_subplot(1, 1, 1)
    for nc in range(len(ncells)):
        for p in range(degrees[0],degrees[-1]+1):

            mask = np.isfinite(timings_i[nc,p-degrees[0]])
            line, = ax.plot(nthreads[mask], timings_i[nc,p-degrees[0]][mask], line_styles[nc],color=colors[p-degrees[0]])

        row = '$n_{{el}}={}^3$'.format(ncells[nc][0][0])
        line, = ax.plot(np.nan*nthreads[mask], np.nan*timings_i[nc,degrees[0]][mask], line_styles[nc],color='k', label=row)

    for p in range(degrees[0],degrees[-1]+1):
        row = '$p={}$'.format(p)
        line, = ax.plot(np.nan*nthreads[mask], np.nan*timings_i[0,degrees[0]][mask],color=colors[p-degrees[0]], label=row)


    box = ax.get_position()
#   ax.set_position([box.x0, box.y0, box.width * 0.3, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_xlabel( xlabel, rotation='horizontal' )
    ax.set_ylabel( r'time [s]' )
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xticks([])
    ax.set_xticks(nthreads)
    ax.set_xticklabels([str(d) for d in nthreads])
    ax.grid(True)
#    ax.title.set_text(title)
    fig.tight_layout()

    fig.savefig("images/"+fname)
