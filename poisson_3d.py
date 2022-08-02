# -*- coding: UTF-8 -*-

from time import time
setup_time1 = time()
import os
import numpy as np

from mpi4py import MPI
from sympy  import pi, sin, symbols

from sympde.calculus import grad, dot
from sympde.topology import ScalarFunctionSpace
from sympde.topology import element_of
from sympde.topology import NormalVector
from sympde.topology import Cube
from sympde.expr     import BilinearForm, LinearForm, integral
from sympde.expr     import Norm
from sympde.expr     import find, EssentialBC

from psydac.api.discretization import discretize
from psydac.api.settings       import PSYDAC_BACKEND_GPYCCEL

x,y,z = symbols('x1, x2, x3')
comm = MPI.COMM_WORLD

def run_poisson_3d(solution, f, ncells, degree, backend=None):

    #+++++++++++++++++++++++++++++++
    # 1. Abstract model
    #+++++++++++++++++++++++++++++++
    domain = Cube()

    B_dirichlet_0 = domain.boundary

    V  = ScalarFunctionSpace('V', domain)
    u  = element_of(V, name='u')
    v  = element_of(V, name='v')

    # Bilinear form a: V x V --> R
    a = BilinearForm((u, v), integral(domain, dot(grad(u), grad(v))))

    # Linear form l: V --> R
    l = LinearForm(v, integral(domain, f * v))

    bc = [EssentialBC(u, 0, domain.boundary)]

    # Variational model
    equation = find(u, forall=v, lhs=a(u, v), rhs=l(v), bc=bc)

    # Error norms
    error  = u - solution
    l2norm = Norm(error, domain, kind='l2')
    h1norm = Norm(error, domain, kind='h1')
    #+++++++++++++++++++++++++++++++
    # 2. Discretization
    #+++++++++++++++++++++++++++++++

    # Create computational domain from topological domain
    domain_h = discretize(domain, ncells=ncells, comm=MPI.COMM_WORLD)

    # Discrete spaces
    Vh = discretize(V, domain_h, degree=degree)

    # Discretize equation using Dirichlet bc
    equation_h = discretize(equation, domain_h, [Vh, Vh], backend=backend)

    # Discretize error norms
    l2norm_h = discretize(l2norm, domain_h, Vh, backend=backend)
    h1norm_h = discretize(h1norm, domain_h, Vh, backend=backend)

    setup_time2 = time()
    T = comm.reduce(setup_time2-setup_time1,op=MPI.MAX)

    infos = {}
    infos['title'] = 'poisson_3d'
    infos['setup_time'] = T
    infos['ncells'] = tuple(ncells)
    infos['degree'] = tuple(degree)
    infos['cart_decomposition'] = tuple(Vh.vector_space.cart.nprocs)
    infos['number_of_threads']  = Vh.vector_space.cart.num_threads

    #+++++++++++++++++++++++++++++++
    # 3. Solution
    #+++++++++++++++++++++++++++++++

    lhs = equation_h.lhs
    rhs = equation_h.rhs
    # Solve linear system
    comm.barrier()
    t1 = time()
    A  = lhs.assemble()
    t2 = time()
    T = comm.reduce(t2-t1,op=MPI.MAX)

    infos['bilinear_form_assembly_time'] = T
    comm.Barrier()
    t1 = time()
    b  = rhs.assemble()
    t2 = time()
    T = comm.reduce(t2-t1,op=MPI.MAX)

    infos['linear_form_assembly_time'] = T

    out = b.copy()
    st  = 0
    for i in range(20):
        comm.Barrier()
        t1 = time()
        A.dot(b, out=out)
        t2 = time()
        b.ghost_regions_in_sync = False
        st += t2-t1

    T = comm.reduce(st/20,op=MPI.MAX)

    infos['dot_product_time'] = T

    st = 0
    for i in range(20):
        comm.Barrier()
        t1 = time()
        b.update_ghost_regions()
        t2 = time()
        st += t2-t1

    T = comm.reduce(st/20,op=MPI.MAX)

    infos['dot_product_communication_time'] = T

    equation_h.set_solver('cg', tol=1e-8, maxiter=3000, info=True)
    equation_h.assemble()
    
    t1 = time()
    # Solve linear system
    uh, info = equation_h.solve()
    t2 = time()

    infos.update(info)
    infos['solve_time'] = comm.reduce(t2-t1,op=MPI.MAX)

    # Compute error norms
    l2_error = l2norm_h.assemble(u=uh)
    h1_error = h1norm_h.assemble(u=uh)

    infos['l2_error'] = l2_error
    infos['h1_error'] = h1_error

    if comm.rank == 0:
        name = (infos['title'],) + infos['ncells'] + infos['degree'] + (comm.size, infos['number_of_threads'])
        name = '_'.join([str(i) for i in name])
        np.save('results/' + name, infos)

#==============================================================================
# 3D Poisson's equation
#==============================================================================
def test_poisson_3d(ncells, degree):

    solution = sin(pi*x)*sin(pi*y)*sin(pi*z)
    f        = 3*pi**2*sin(pi*x)*sin(pi*y)*sin(pi*z)

    run_poisson_3d(solution, f, ncells=ncells, degree=degree, backend=PSYDAC_BACKEND_GPYCCEL)


#==============================================================================
# Parser
#==============================================================================
def parse_input_arguments():

    import argparse

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description     = "Solve Poisson's equation on a 3D domain."
    )

    parser.add_argument( '-d',
        type    = int,
        nargs   = 3,
        default = [2,2,2],
        metavar = ('P1','P2','P3'),
        dest    = 'degree',
        help    = 'Spline degree along each dimension'
    )

    parser.add_argument( '-n',
        type    = int,
        nargs   = 3,
        default = [10,10,10],
        metavar = ('N1','N2','N3'),
        dest    = 'ncells',
        help    = 'Number of grid cells (elements) along each dimension'
    )


    return parser.parse_args()

#==============================================================================
# Script functionality
#==============================================================================
if __name__ == '__main__':

    args = parse_input_arguments()
    test_poisson_3d(**vars(args))


