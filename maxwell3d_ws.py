# -*- coding: UTF-8 -*-

from time import time
setup_time1 = time()

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
from psydac.linalg.utilities         import array_to_psydac

PSYDAC_BACKEND_GPYCCEL = PSYDAC_BACKEND_GPYCCEL.copy()

if int(os.environ.get('OMP_NUM_THREADS', 1))>1:
    PSYDAC_BACKEND_GPYCCEL['openmp'] = True


def remove_folder(path):
    os.system('rm -rf "%s" &' % path)
#===============================================================================
def union(name, domains):
    assert len(domains)>1
    domain = domains[0]
    for p in domains[1:]:
        domain = domain.join(p, name=name)
    return domain

def set_interfaces(domain, interfaces):
    for I in interfaces:
        domain = domain.join(domain, domain.name, bnd_minus=I[0], bnd_plus=I[1])
    return domain

def join(name, domains, interfaces):
    domain = union(name, domains)
    domain = set_interfaces(domain, interfaces)
    return domain

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
    a2 = BilinearForm((E,TE), integral(domain, dot(E, curl(TE)))- integral(boundary, expr_b) + integral(I, expr_I))
    a3 = BilinearForm((E,TE), -integral(domain, dot(E, curl(TE)))- integral(I, expr_I))

    l1 = LinearForm(TE, integral(domain, dot(Eex, TE)))
    l2 = LinearForm(TE, integral(domain, dot(Bex, TE)))
    l3 = LinearForm(TE, integral(domain, dot(J, TE)))
    l4 = LinearForm(TE, integral(boundary, dot(cross(Eex, nn),TE)))

    l2norm = Norm(error, domain, kind='l2')
    #+++++++++++++++++++++++++++++++
    # 2. Discretization
    #+++++++++++++++++++++++++++++++

    domain_h = discretize(domain, ncells=ncells, comm=comm)
    Vh       = discretize(V, domain_h, degree=degree)

    a1_h = discretize(a1, domain_h, [Vh, Vh], backend=backend)

    l1_h = discretize(l1, domain_h, Vh, backend=backend)

    comm.Barrier()

    setup_time2 = time()
    tt = comm.reduce(setup_time2-setup_time1,op=MPI.MAX)

    # store setup info
    infos = {}
    infos['title'] = 'maxwell_3d'
    infos['setup_time'] = tt
    infos['ncells'] = tuple(ncells)
    infos['degree'] = tuple(degree)
    infos['cg_niter'] = cg_niter
    infos['cart_decompositions'] = tuple(tuple(i.nprocs) for i in domain_h.ddm.domains)
    infos['number_of_threads']  = domain_h.ddm.num_threads

    comm.Barrier()

    t1 = time()
    M1 = a1_h.assemble()
    t2 = time()

    tt = comm.reduce(t2-t1,op=MPI.MAX)

    # store bilinear forms assembly timing
    infos['bilinear_form_assembly_time'] = tt

    b   = l1_h.assemble(t=0.)
    out = b.copy()
    M1.dot(b)
    b.ghost_regions_in_sync = False
    st  = 0
    for i in range(200):
        comm.Barrier()
        t1 = time()
        M1.dot(b, out=out)
        t2 = time()
        b.ghost_regions_in_sync = False
        st += t2-t1

    tt = comm.reduce(st/200,op=MPI.MAX)

    # store mass matrix dot product timing
    infos['dot_product_time'] = tt

    st = 0
    for i in range(20):
        comm.Barrier()
        t1 = time()
        b.update_ghost_regions()
        t2 = time()
        st += t2-t1

    tt = comm.reduce(st/20,op=MPI.MAX)

    # store communication timing
    infos['dot_product_communication_time'] = tt

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


    A1 = Cube('A1',bounds1=(0, 0.5), bounds2=(0., 0.5), bounds3=(0., 0.5))
    A2 = Cube('A2',bounds1=(0.5, 1), bounds2=(0., 0.5), bounds3=(0.,0.5))
    A3 = Cube('A3',bounds1=(0, 0.5), bounds2=(0.5, 1), bounds3=(0, 0.5))
    A4 = Cube('A4',bounds1=(0.5, 1.), bounds2=(0.5, 1), bounds3=(0, 0.5))
    A5 = Cube('A5',bounds1=(0, 0.5), bounds2=(0, 0.5), bounds3=(0.5, 1))
    A6 = Cube('A6',bounds1=(0.5, 1), bounds2=(0, 0.5), bounds3=(0.5, 1))
    A7 = Cube('A7',bounds1=(0, 0.5), bounds2=(0.5, 1), bounds3=(0.5, 1))

    patches = [A1, A2, A3, A4, A5, A6, A7]
    interfaces = [
                 [A1.get_boundary(axis=0, ext=1), A2.get_boundary(axis=0,ext=-1)],
                 [A3.get_boundary(axis=0, ext=1), A4.get_boundary(axis=0,ext=-1)],
                 [A5.get_boundary(axis=0, ext=1), A6.get_boundary(axis=0,ext=-1)],

                 [A1.get_boundary(axis=1, ext=1), A3.get_boundary(axis=1,ext=-1)],
                 [A2.get_boundary(axis=1, ext=1), A4.get_boundary(axis=1,ext=-1)],
                 [A5.get_boundary(axis=1, ext=1), A7.get_boundary(axis=1,ext=-1)],

                 [A1.get_boundary(axis=2, ext=1), A5.get_boundary(axis=2,ext=-1)],
                 [A2.get_boundary(axis=2, ext=1), A6.get_boundary(axis=2,ext=-1)],
                 [A3.get_boundary(axis=2, ext=1), A7.get_boundary(axis=2,ext=-1)],
                 ]

    domain = join('Omega', patches, interfaces)

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

