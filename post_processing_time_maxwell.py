import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})

results_folder = 'results/'
number_of_mpi_procs = [1*35,2*35,4*35,8*35,16*35,32*35,64*35,128*35]
number_of_threads   = 1
filename            = 'maxwell3d'

number_of_mpi_procs = [32*7,32*7*2,32*7*2**2,32*7*2**3,32*7*2**4]
ncells          = [80]
degrees         = [3,4,5]
number_of_threads = 1

timmings_bi_assembly     = np.zeros((len(ncells),len(degrees), len(number_of_mpi_procs)))
timmings_time_integrator = np.zeros((len(ncells),len(degrees), len(number_of_mpi_procs)))
timmings_dot_p           = np.zeros((len(ncells),len(degrees), len(number_of_mpi_procs)))
scaling_bi_assembly      = timmings_bi_assembly.copy()
scaling_time_integrator  = timmings_time_integrator.copy()
scaling_dot_p            = timmings_dot_p.copy()

for i1,nc in enumerate(ncells):
    for i2,d in enumerate(degrees):
        for i3, mpi_p in enumerate(number_of_mpi_procs):
            names = ('maxwell_3d',) + (nc,)*3 +(d,)*3 + (mpi_p, number_of_threads)
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

        k = [j for j in range(len(number_of_mpi_procs)) if not np.isnan(timmings_bi_assembly[i1,i2,j])] + [0]
        nn = [np.nan]*k[0] + [2**(i-k[0]) for i in range(k[0],len(number_of_mpi_procs))]
        scaling_bi_assembly[i1,i2,:]     = timmings_bi_assembly[i1,i2,k[0]]/timmings_bi_assembly[i1,i2,:]/nn
        scaling_time_integrator[i1,i2,:] = timmings_time_integrator[i1,i2,k[0]]/timmings_time_integrator[i1,i2,:]/nn
        scaling_dot_p[i1,i2,:]           = timmings_dot_p[i1,i2,k[0]]/timmings_dot_p[i1,i2,:]/nn

number_of_mpi_procs = [2*7,2*7*2,2*7*2**2,2*7*2**3,2*7*2**4]
number_of_threads = 16
timmings_bi_assembly_mth     = np.zeros((len(ncells),len(degrees), len(number_of_mpi_procs)))
timmings_time_integrator_mth = np.zeros((len(ncells),len(degrees), len(number_of_mpi_procs)))
timmings_dot_p_mth           = np.zeros((len(ncells),len(degrees), len(number_of_mpi_procs)))
scaling_bi_assembly_mth      = timmings_bi_assembly_mth.copy()
scaling_time_integrator_mth  = timmings_time_integrator_mth.copy()
scaling_dot_p_mth            = timmings_dot_p_mth.copy()
for i1,nc in enumerate(ncells):
    for i2,d in enumerate(degrees):
        for i3, mpi_p in enumerate(number_of_mpi_procs):
            names = ('maxwell_3d',) + (nc,)*3 +(d,)*3 + (mpi_p, number_of_threads)
            filename = '_'.join([str(i) for i in names])
            try:
                T = np.load(results_folder+filename+'.npy', allow_pickle=True)
                T = T.item()
                timmings_bi_assembly_mth[i1,i2,i3]     = min(T['bilinear_form_assembly_time'],T.get('bilinear_form_assembly_time2', T['bilinear_form_assembly_time']))
#                timmings_time_integrator_mth[i1,i2,i3] = T['solution']
                timmings_dot_p_mth[i1,i2,i3]           = T['dot_product_time']
            except:
                timmings_bi_assembly_mth[i1,i2,i3] = np.nan
                timmings_time_integrator_mth[i1,i2,i3] = np.nan
                timmings_dot_p_mth[i1,i2,i3] = np.nan

        k = [j for j in range(len(number_of_mpi_procs)) if not np.isnan(timmings_bi_assembly_mth[i1,i2,j])] + [0]
        nn = [np.nan]*k[0] + [2**(i-k[0]) for i in range(k[0],len(number_of_mpi_procs))]
        scaling_bi_assembly_mth[i1,i2,:]     = timmings_bi_assembly_mth[i1,i2,k[0]]/timmings_bi_assembly_mth[i1,i2,:]/nn
        scaling_time_integrator_mth[i1,i2,:] = timmings_time_integrator_mth[i1,i2,k[0]]/timmings_time_integrator_mth[i1,i2,:]/nn
        scaling_dot_p_mth[i1,i2,:]           = timmings_dot_p_mth[i1,i2,k[0]]/timmings_dot_p_mth[i1,i2,:]/nn

number_of_mpi_procs = np.array([32*7,32*7*2,32*7*2**2,32*7*2**3,32*7*2**4])
##-----------------------------------------------------------------------------------------------------------------------
#from tabulate import tabulate
#headers = [""] + [str(np) for np in number_of_mpi_procs]

#if not all(np.isnan(v) for v in timmings_bi_assembly.flatten()):
#    print("="*45,"Timings of the Matrix Assembly of Time dependent Maxwell", "="*45)
#    T = np.around(scaling_bi_assembly_mth, decimals=5)
#    newT = []
#    for i1,nc in enumerate(ncells):
#        for i2,d in enumerate(degrees):
#            newT.append(["nc = {} ** 3 , p = {}".format(nc,d)] +  T[i1,i2].tolist())
#        newT.append(["   "]*len(T[0]))

#    print(tabulate(newT, headers=headers, tablefmt="grid"))
#    print("\n")
#raise

#====================================================================================================
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
fnames = ['matrix_assembly_time_maxwell_strong_scaling', 'matrix_vector_product_time_maxwell_strong_scaling','matrix_assembly_time_maxwell_strong_scaling_multi_threading', 'matrix_vector_product_time_maxwell_strong_scaling_multi_threading']
xaxist = [r'Number of nodes']*4
timings = [[timmings_bi_assembly,timmings_bi_assembly_mth], [timmings_dot_p, timmings_dot_p_mth]]
number_of_nodes = np.array([7*1,7*2,7*4,7*8,7*16])

for title,fname,timings_i,xlabel in zip(titles, fnames, timings,xaxist):
    fig = plt.figure(figsize=(15,8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(number_of_nodes,[5*np.nanmax(timings_i)/2**d for d in range(len(number_of_nodes))], color='black', linestyle='dashed', label='perfect scaling')
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
    ax.set_xlabel(xlabel, rotation='horizontal' )
    ax.set_ylabel( r'time [s]' )
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xticks([])
    ax.set_xticks(number_of_nodes)
    ax.set_xticklabels([str(d) for d in number_of_nodes])
    ax.grid(True)
#    ax.title.set_text(title)
    fig.tight_layout(rect=[0, 0.05, 1, 1])

    fig.savefig("images/"+fname)
