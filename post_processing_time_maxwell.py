import numpy as np

results_folder = 'results/'
ncells   = [125,150,175]
degrees  = [2,3]
number_of_mpi_procs = [1*35,2*35,4*35,8*35,16*35,32*35,64*35,128*35]
number_of_threads   = 1
filename            = 'maxwell3d'

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

#-----------------------------------------------------------------------------------------------------------------------
from tabulate import tabulate
headers = [""] + [str(np) for np in number_of_mpi_procs]

if not all(np.isnan(v) for v in timmings_bi_assembly.flatten()):
    print("="*45,"Timings of the Matrix Assembly of Time dependent Maxwell", "="*45)
    T = np.around(scaling_bi_assembly, decimals=5)
    newT = []
    for i1,nc in enumerate(ncells):
        for i2,d in enumerate(degrees):
            newT.append(["nc = {} ** 3 , p = {}".format(nc,d)] +  T[i1,i2].tolist())
        newT.append(["   "]*len(T[0]))

    print(tabulate(newT, headers=headers, tablefmt="grid"))
    print("\n")

if not all(np.isnan(v) for v in timmings_time_integrator.flatten()):
    print("="*45,"Timings of the time intergrator of Time dependent Maxwell", "="*45)
    T = np.around(scaling_time_integrator, decimals=5)
    newT = []
    for i1,nc in enumerate(ncells):
        for i2,d in enumerate(degrees):
            newT.append(["nc = {} ** 3 , p = {}".format(nc,d)] +  T[i1,i2].tolist())
        newT.append(["   "]*len(T[0]))

    print(tabulate(newT, headers=headers, tablefmt="grid"))
    print("\n")

if not all(np.isnan(v) for v in timmings_dot_p.flatten()):
    print("="*45,"Timings of the Matrix vector product of the Hcurl Mass Matrix", "="*45)
    T = np.around(scaling_dot_p, decimals=5)
    newT = []
    for i1,nc in enumerate(ncells):
        for i2,d in enumerate(degrees):
            newT.append(["nc = {} ** 3 , p = {}".format(nc,d)] +  T[i1,i2].tolist())
        newT.append(["   "]*len(T[0]))

    print(tabulate(newT, headers=headers, tablefmt="grid"))
    print("\n")

#====================================================================================================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.legend_handler import HandlerLine2D

from itertools import product

colors = cm.rainbow(np.linspace(0, 1, len(degrees)*len(ncells)).reshape((len(ncells),len(degrees))))
titles = ['Matrix Assembly']
names  = ['assembly']
timings = timmings_bi_assembly
for title,name in zip(titles, names):
    fig = plt.figure()
#    fig.suptitle(title)
    for nc in range(len(ncells)):
        ax = fig.add_subplot(2, 2, nc+1)
        for p in range(2,4):
            row = '$p={}$'.format(p)
            line, = ax.plot(number_of_mpi_procs, timings[nc,p-2], 's-',color=colors[nc, p-2], label='Psydac,' + row)
            #ax.plot(number_of_mpi_procs,timings[nc,p-2],'s', color=colors[nc,p-2])

        ax.plot(number_of_mpi_procs,[0.3*np.nanmax(timings[0,p-2])/2**d for d in range(len(number_of_mpi_procs))], color='black', linestyle='dashed', label='perfect scaling')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        ax.legend(handler_map={line: HandlerLine2D(numpoints=4)})
        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        #ax.legend(bbox_to_anchor=(1, 1), loc='upper right', ncol=1, fontsize=7)
        ax.set_xlabel( r'nprocs', rotation='horizontal' )
        ax.set_ylabel( r'time [s]' )
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xticks([])
        ax.set_xticks(number_of_mpi_procs)
        ax.set_xticklabels([str(d) for d in number_of_mpi_procs])
        ax.grid(True)
        ax.title.set_text('$ncells={}^3$'.format(ncells[nc]))
    fig.subplots_adjust(hspace = .4, wspace=.4, left=0.05, right=0.9)
    #fig.tight_layout()
plt.show()
