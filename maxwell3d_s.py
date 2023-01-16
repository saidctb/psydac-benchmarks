# -*- coding: UTF-8 -*-

from time import time

import os
import shutil
import numpy as np
from mpi4py import MPI
from sympy  import pi, sin, cos, Tuple, ImmutableDenseMatrix as Matrix

from scipy.sparse.linalg import spsolve, inv

from sympde.calculus        import grad, dot, curl, cross
from sympde.calculus        import minus, plus
from sympde.topology        import VectorFunctionSpace
from sympde.topology        import elements_of
from sympde.topology        import NormalVector
from sympde.topology        import Square, Cube
from sympde.topology        import IdentityMapping, PolarMapping
from sympde.expr.expr       import LinearForm, BilinearForm
from sympde.expr.expr       import integral
from sympde.expr            import TerminalExpr
from sympde.expr.expr       import Norm
from sympde.expr.equation   import find, EssentialBC
from sympde.utilities.utils import plot_domain
from sympde.core            import Constant

from psydac.api.discretization       import discretize
from psydac.api.tests.build_domain   import build_pretzel
from psydac.fem.basic                import FemField
from psydac.linalg.iterative_solvers import *
from psydac.api.settings             import PSYDAC_BACKEND_GPYCCEL,PSYDAC_DEFAULT_FOLDER
from psydac.feec.pull_push           import pull_2d_hcurl
from psydac.fem.vector               import ProductFemSpace
from psydac.linalg.iterative_solvers import pcg
from psydac.api.postprocessing       import OutputManager
from psydac.linalg.utilities         import array_to_psydac

PSYDAC_BACKEND_GPYCCEL = PSYDAC_BACKEND_GPYCCEL.copy()

if int(os.environ.get('OMP_NUM_THREADS', 1))>1:
    PSYDAC_BACKEND_GPYCCEL['openmp'] = True

from hcurl_mass_matrix import assemble_matrix
#==============================================================================
def run_maxwell_3d(Eex, Bex, J, domain, ncells, degree,  dt, niter, T, backend, comm, cg_niter,store=None):

    backend['folder'] = "maxwell_3d_psydac_{}_{}_{}_{}".format(ncells[0], degree[0], comm.size, int(os.environ.get('OMP_NUM_THREADS', 1)))
    backend['flags']  = "-O3 -march=native -mtune=native  -mavx -ffast-math -ffree-line-length-none"
    PSYDAC_DEFAULT_FOLDER['name'] = '__psydac__' + backend['folder']
    #+++++++++++++++++++++++++++++++
    # 1. Abstract model
    #+++++++++++++++++++++++++++++++
    V  = VectorFunctionSpace('V', domain, kind='hcurl')

    E, TE  = elements_of(V, names='E, TE')
    nn = NormalVector('nn')

    error    = Matrix([E[0]-Eex[0],E[1]-Eex[1],E[2]-Eex[2]])
    I        = domain.interfaces
    boundary = domain.boundary

    jump  = lambda w:minus(w)-plus(w)
    avr   = lambda w:0.5*plus(w) + 0.5*minus(w)

    expr_I  =  dot(cross(jump(TE),nn),avr(E))
    expr_b  =  dot(cross(E, nn),TE)

    # Bilinear form a: V x V --> R
    a1 = BilinearForm((E,TE), integral(domain, dot(E, TE)))

    domain_h = discretize(domain, ncells=ncells, comm=comm)
    Vh       = discretize(V, domain_h, degree=degree)

    a1_h = discretize(a1, domain_h, [Vh, Vh], backend=backend)

    # store setup info
    infos = {}
    infos['title'] = 'maxwell_3d_s'
    infos['ncells'] = tuple(ncells)
    infos['degree'] = tuple(degree)
    infos['cg_niter'] = cg_niter
    infos['cart_decompositions'] = (32,)
    infos['number_of_threads']  = domain_h.ddm.num_threads
    if infos['number_of_threads']>1:
        a1_h._func = assemble_matrix
    comm.Barrier()

    t1 = time()
    M1 = a1_h.assemble()
    t2 = time()

    tt = comm.reduce(t2-t1,op=MPI.MAX)

    # store bilinear forms assembly timing
    infos['bilinear_form_assembly_time'] = tt

    comm.Barrier()

    t1 = time()
    M1 = a1_h.assemble()
    t2 = time()

    tt = comm.reduce(t2-t1,op=MPI.MAX)

    # store bilinear forms assembly timing
    infos['bilinear_form_assembly_time2'] = tt


    if comm.rank == 0:
        name = (infos['title'],) + infos['ncells'] + infos['degree'] + (comm.size, infos['number_of_threads'])
        name = '_'.join([str(i) for i in name])
        np.save('results/' + name, infos)


    return infos

def parse_input_arguments():

    import argparse

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description     = "Solve the Time dependent Maxwell equation on a 3D domain with" +
                          " inhomogeneous Dirichlet boundary conditions."
    )

    parser.add_argument( '-d',
        type    = int,
        nargs   = 3,
        default = [2, 2, 2],
        metavar = ('P1', 'P2', 'P3'),
        dest    = 'degree',
        help    = 'Spline degree along each dimension'
    )

    parser.add_argument( '-n',
        type    = int,
        nargs   = 3,
        default = [10, 10, 10],
        metavar = ('N1', 'N2', 'N3'),
        dest    = 'ncells',
        help    = 'Number of grid cells (elements) along each dimension'
    )

#    parser.add_argument( '-ns',
#        type    = int,
#        default = 10,
#        dest    = 'NSTEPS',
#        metavar = 'NSTEPS',
#        help    = 'Number of time-steps to be taken'
#    )

#    parser.add_argument( '-T',
#        type    = float,
#        dest    = 'END_TIME',
#        metavar = 'END_TIME',
#        help    = 'Run simulation until given final time'
#    )
#    # ...

    parser.add_argument( '-s',
        action  = 'store_true',
        dest    = 'store',
        help    = 'Save output files'
    )

    parser.add_argument('--a', action='store_true', \
                       help='Use analytical mapping.', \
                       dest='analytical')

    parser.add_argument( '-m',
        type    = str,
        nargs   = 1,
        default = ['identity'],
        dest    = 'mapping',
        help    = 'mapping'
    )
    return parser.parse_args()
    

if __name__ == '__main__':

    from psydac.api.postprocessing import PostProcessManager

    domain = Cube('A1',bounds1=(0, 1), bounds2=(0., 1), bounds3=(0., 1))

    x,y,z = domain.coordinates
    t     = Constant(name='t')

    Eex   = Matrix([sin(pi*y)*sin(pi*z)*cos(t),
                  sin(pi*x)*sin(pi*z)*cos(t),
                  sin(pi*x)*sin(pi*y)*cos(t)])

    Bex   = Matrix([
            [pi*(-cos(pi*y) + cos(pi*z))*sin(pi*x)*sin(t)],
            [pi*( cos(pi*x) - cos(pi*z))*sin(pi*y)*sin(t)],
            [pi*(-cos(pi*x) + cos(pi*y))*sin(pi*z)*sin(t)]])

    J = Eex.diff(t) - TerminalExpr(curl(Bex), domain)

    args = vars(parse_input_arguments())

    ncells = args['ncells']
    degree = args['degree']

    #time parameters
    niter = 10
    T     = 0.01
    dt    = T/niter
    cg_niter = 10
 
    comm  = MPI.COMM_WORLD
    store = args['store']

    backend = PSYDAC_BACKEND_GPYCCEL

    run_maxwell_3d(Eex, Bex, J, domain, ncells, degree,  dt, niter, T, backend, comm, cg_niter,store=store)

