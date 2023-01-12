import numpy as np

results_folder = 'results/'
filename       = 'maxwell3d'

#ncells0 = [[10,10,10],[20,10,10],[20,20,10],[20,20,20],[40,20,20],[40,40,20],[40,40,40],[80,40,40],[80,80,40],[80,80,80]]
ncells1 = [[11,11,11],[22,11,11],[22,22,11],[22,22,22],[44,22,22],[44,44,22],[44,44,44],[88,44,44],[88,88,44],[88,88,88]]
ncells2 = [[12,12,12],[24,12,12],[24,24,12],[24,24,24],[48,24,24],[48,48,24],[48,48,48],[96,48,48],[96,96,48],[96,96,96]]
ncells = [ncells1,ncells2]
degrees = [3,4,5]

nnodes          = [1,1,1,2,4,7,7*2,7*2**2,7*2**3,7*2**4]
ntasks_per_node = [7,14,28,28,28,32,32,32,32,32]
nthreads        = [1,1,1,1,1,1,1,1,1,1,1,1]

timmings_bi_assembly     = np.zeros((len(ncells),len(degrees), len(nnodes)))
timmings_time_integrator = np.zeros((len(ncells),len(degrees), len(nnodes)))
timmings_dot_p           = np.zeros((len(ncells),len(degrees), len(nnodes)))

for i1,nc in enumerate(ncells):
    for i2,d in enumerate(degrees):
        for i3, (nn,nt,nth) in enumerate(zip(nnodes, ntasks_per_node, nthreads)):
            mpi_p = nn*nt
            names = ('maxwell_3d',) + (nc[i3][0],nc[i3][1],nc[i3][2]) +(d,)*3 + (mpi_p, nth)
            filename = '_'.join([str(i) for i in names])
            try:
                T = np.load(results_folder+filename+'.npy', allow_pickle=True)
                T = T.item()
                timmings_bi_assembly[i1,i2,i3]     = min(T['bilinear_form_assembly_time'],T.get('blinear_form_assembly_time2', T['bilinear_form_assembly_time']))
                timmings_time_integrator[i1,i2,i3] = T['solution']
                timmings_dot_p[i1,i2,i3]           = T['dot_product_time']
            except:
                timmings_bi_assembly[i1,i2,i3] = np.nan
                timmings_time_integrator[i1,i2,i3] = np.nan
                timmings_dot_p[i1,i2,i3] = np.nan

nnodes          = [7,7,7,7,7,7,7*2,7*2**2,7*2**3,7*2**4]
ntasks_per_node = [1,1,1,1,1,2,2,2,2,2]
nthreads        = [1,2,4,8,16,16,16,16,16,16]

timmings_bi_assembly_mth     = np.zeros((len(ncells),len(degrees), len(nnodes)))
timmings_time_integrator_mth = np.zeros((len(ncells),len(degrees), len(nnodes)))
timmings_dot_p_mth           = np.zeros((len(ncells),len(degrees), len(nnodes)))

for i1,nc in enumerate(ncells):
    for i2,d in enumerate(degrees):
        for i3, (nn,nt,nth) in enumerate(zip(nnodes, ntasks_per_node, nthreads)):
            mpi_p = nn*nt
            names = ('maxwell_3d',) + (nc[i3][0],nc[i3][1],nc[i3][2]) +(d,)*3 + (mpi_p, nth)
            filename = '_'.join([str(i) for i in names])

            try:
                T = np.load(results_folder+filename+'.npy', allow_pickle=True)
                T = T.item()
                timmings_bi_assembly_mth[i1,i2,i3]     = min(T['bilinear_form_assembly_time'],T.get('blinear_form_assembly_time2', T['bilinear_form_assembly_time']))
                timmings_time_integrator_mth[i1,i2,i3] = T['solution']
                timmings_dot_p_mth[i1,i2,i3]           = T['dot_product_time']
            except:
                timmings_bi_assembly_mth[i1,i2,i3] = np.nan
                timmings_time_integrator_mth[i1,i2,i3] = np.nan
                timmings_dot_p_mth[i1,i2,i3] = np.nan
##-----------------------------------------------------------------------------------------------------------------------
#from tabulate import tabulate
#headers = [""] + [str(nn*nt) for nn,nt in zip(nnodes, ntasks_per_node)]

#if not all(np.isnan(v) for v in timmings_bi_assembly.flatten()):
#    print("="*45,"Timings of the Matrix Assembly of Time dependent Maxwell", "="*45)
#    T = np.around(timmings_bi_assembly, decimals=5)
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
fnames = ['matrix_assembly_tm_ws', 'matrix_vector_product_tm_ws','matrix_assembly_mth_tm_ws', 'matrix_vector_product_mth_tm_ws']
xaxist = [r'number of mpi procs', r'number of mpi procs',r'number of threads',r'number of threads']
timings = [timmings_bi_assembly, timmings_dot_p, timmings_bi_assembly_mth, timmings_dot_p_mth]

nnodes   = [1,1,1,2,4,7,7*2,7*2**2,7*2**3,7*2**4]
nthreads = np.array([nn*nt*nth for nn,nt,nth in zip(nnodes, ntasks_per_node, nthreads)])

for title,fname,timings_i,xlabel in zip(titles, fnames, timings,xaxist):
    fig = plt.figure(figsize=(10,15))
    ax = fig.add_subplot(1, 1, 1)
    for nc in range(len(ncells)):
        for p in range(degrees[0],degrees[-1]+1):

            mask = np.isfinite(timings_i[nc,p-degrees[0]])
            line, = ax.plot(nthreads[mask], timings_i[nc,p-degrees[0]][mask], line_styles[nc],color=colors[p-degrees[0]])

        row = '$n_{{el}}={}^3$'.format(ncells[nc][0][0])
        line, = ax.plot(np.nan*nthreads[mask], np.nan*timings_i[nc,p-degrees[0]][mask], line_styles[nc],color='k', label=row)

    for p in range(degrees[0],degrees[-1]+1):
        row = '$p={}$'.format(p)
        line, = ax.plot(np.nan*nthreads[mask], np.nan*timings_i[0,p-degrees[0]][mask],color=colors[p-degrees[0]], label=row)


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

    fig.savefig(fname)
