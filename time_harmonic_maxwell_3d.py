# -*- coding: UTF-8 -*-
from time import time
setup_time1 = time()
import os
import shutil
import numpy as np

from mpi4py import MPI
from sympy import pi, cos, sin, sqrt, Matrix, Tuple, symbols

from sympde.calculus import grad, dot, inner, div, curl, cross
from sympde.topology import NormalVector
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import ProductSpace
from sympde.topology import element_of, elements_of
from sympde.topology import Cube
from sympde.expr     import BilinearForm, LinearForm, integral
from sympde.expr     import Norm
from sympde.expr     import find, EssentialBC

from psydac.api.discretization import discretize
from psydac.api.settings       import PSYDAC_BACKEND_GPYCCEL

x,y,z    = symbols('x1, x2, x3')
comm = MPI.COMM_WORLD

def remove_folder(path):
    shutil.rmtree(path)

def run_maxwell_time_harmonic_3d(uex, f, alpha, ncells, degree, backend):

    backend['folder'] = "time_harmonic_3d_psydac_{}_{}_{}".format(ncells[0], degree[0], comm.size)
    backend['flags']  = "-O3 -march=native -mtune=native  -mavx -ffast-math"

    # ... abstract model
    domain = Cube('A')

    V  = VectorFunctionSpace('V', domain, kind='hcurl')

    u  = element_of(V, name='u')
    v  = element_of(V, name='v')
    nn   = NormalVector('nn')

    # Bilinear form a: V x V --> R
    a   = BilinearForm((u, v), integral(domain, dot(curl(u),curl(v)) + alpha*dot(u,v)) + integral(domain.boundary, 1e10 * dot(cross(u, nn),cross(v, nn))) )

    # Linear form l: V --> R
    l = LinearForm(v, integral(domain, dot(f,v)))

    # Variational model
    equation = find(u, forall=v, lhs=a(u, v), rhs=l(v))

    # l2 error
    error   = Matrix([u[0]-uex[0],u[1]-uex[1], u[2]-uex[2]])
    l2norm  = Norm(error, domain, kind='l2')

    #+++++++++++++++++++++++++++++++
    # 2. Discretization
    #+++++++++++++++++++++++++++++++

    # Create computational domain from topological domain
    domain_h = discretize(domain, ncells=ncells, comm=comm)

    # Discrete spaces
    Vh = discretize(V, domain_h, degree=degree)

    # Discretize equation
    equation_h  = discretize(equation, domain_h, [Vh, Vh], backend=backend)

    l2_norm_h = discretize(l2norm, domain_h, Vh, backend=backend)

    if comm.rank == 0:
        try:
            remove_folder(backend['folder'])
        except:
            pass
    comm.Barrier()

    setup_time2 = time()
    T = comm.reduce(setup_time2-setup_time1,op=MPI.MAX)

    infos = {}
    infos['title'] = 'time_harmonic_maxwell_3d'
    infos['setup_time'] = T
    infos['ncells'] = tuple(ncells)
    infos['degree'] = tuple(degree)
    infos['cart_decomposition'] = tuple(Vh.spaces[0].vector_space.cart.nprocs)
    infos['number_of_threads']  = Vh.spaces[0].vector_space.cart.num_threads

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
 
    equation_h.set_solver('pcg', pc='jacobi', tol=1e-8, maxiter=3000, info=True)
    equation_h.assemble()

    # Solve linear system
    t1 = time()
    uh, info = equation_h.solve()
    t2 = time()

    infos.update(info)
    infos['solve_time'] = comm.reduce(t2-t1,op=MPI.MAX)


    # Compute error norms
    l2_error = l2_norm_h.assemble(u=uh)

    infos['l2_error'] = l2_error
    if comm.rank == 0:
        name = (infos['title'],) + infos['ncells'] + infos['degree'] + (comm.size, infos['number_of_threads'])
        name = '_'.join([str(i) for i in name])
        np.save('results/' + name, infos)

def test_maxwell_time_harmonic_3d(ncells, degree):

    alpha    = 1.
    uex      = uex = Tuple(sin(pi*y)*sin(pi*z), sin(pi*x)*cos(pi*y)*sin(pi*z), sin(pi*x)*cos(pi*z)*sin(pi*y))
    f        = Tuple(alpha*sin(pi*y)*sin(pi*z) -2*pi**2*sin(pi*y)*sin(pi*z)*cos(pi*x) + 2*pi**2*sin(pi*y)*sin(pi*z),
                     alpha*sin(pi*x)*cos(pi*y)*sin(pi*z) + pi**2*sin(pi*x)*sin(pi*z)*cos(pi*y),
                     alpha*sin(pi*x)*cos(pi*z)*sin(pi*y) + pi**2*sin(pi*x)*sin(pi*y)*cos(pi*z))

    run_maxwell_time_harmonic_3d(uex, f, alpha, ncells=ncells, degree=degree, backend=PSYDAC_BACKEND_GPYCCEL)

#==============================================================================
# Parser
#==============================================================================
def parse_input_arguments():

    import argparse

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description     = "Solve Time harmonic maxwell's equation on a 3D domain."
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
    test_maxwell_time_harmonic_3d(**vars(args))
