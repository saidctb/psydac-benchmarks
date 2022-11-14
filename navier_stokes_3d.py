# -*- coding: UTF-8 -*-

import os
import pytest
import numpy as np
from sympy import pi, cos, sin, sqrt, exp, ImmutableDenseMatrix as Matrix, Tuple, lambdify
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import gmres as sp_gmres
from scipy.sparse.linalg import minres as sp_minres
from scipy.sparse.linalg import cg as sp_cg
from scipy.sparse.linalg import bicg as sp_bicg
from scipy.sparse.linalg import bicgstab as sp_bicgstab

from sympde.calculus import grad, dot, inner, div, curl, cross
from sympde.calculus import Transpose, laplace
from sympde.topology import NormalVector, Mapping, Cube
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import ProductSpace
from sympde.topology import element_of, elements_of
from sympde.topology import Domain, Square, Union
from sympde.expr     import BilinearForm, LinearForm, integral
from sympde.expr     import Norm
from sympde.expr     import find, EssentialBC
from sympde.core     import Constant
from sympde.expr     import TerminalExpr
from sympde.expr     import linearize

from psydac.api.essential_bc   import apply_essential_bc
from psydac.api.postprocessing import OutputManager
from psydac.fem.basic          import FemField
from psydac.fem.vector         import ProductFemSpace
from psydac.api.discretization import discretize
from psydac.linalg.utilities   import petsc_to_psydac
from psydac.linalg.stencil     import *
from psydac.linalg.block       import *
from psydac.api.settings       import PSYDAC_BACKEND_GPYCCEL
from psydac.utilities.utils    import refine_array_1d, animate_field, split_space, split_field
from psydac.linalg.iterative_solvers import cg, pcg, bicg, lsmr

from petsc4py import PETSc
from mpi4py import MPI
comm = MPI.COMM_WORLD

ksp = PETSc.KSP()
ksp.create(comm)
tol = 1e-11

opts = PETSc.Options()
opts["ksp_type"] = "gmres"
opts["pc_type"] = "none"
opts["ksp_rtol"] = tol
opts["ksp_atol"] = tol
ksp.setFromOptions()
#==============================================================================
class ConcentricPipe(Mapping): 
    """
    Represents a Concentric Pipe Mapping object.
    """ 
    _expressions = {'x': 'x1', 
                    'y': 'c1 + (rmin*(1-x2)+rmax*x2)*cos(x3)', 
                    'z': 'c2 + (rmin*(1-x2)+rmax*x2)*sin(x3)'} 
 
    _ldim        = 3 
    _pdim        = 3

#==============================================================================

#------------------------------------------------------------------------------
def petsc_solver(M, b):
    space = M.domain

    M = M.topetsc()
    b = b.topetsc()

    ksp.setOperators(M)
    x = b.duplicate()
    ksp.solve(b, x)
    x  = petsc_to_psydac(x, space)
    return x,0

#------------------------------------------------------------------------------
def psydac_solver(M, b):
    return lsmr(M, M.T, b, maxiter=10000, tol=1e-6)

#==============================================================================
def run_time_dependent_navier_stokes_2d(domain, ncells, degree, multiplicity, boundary_h, boundary_i, dt_h, nt, newton_tol=1e-4, max_newton_iter=100, petsc=False, store=None):

    # ... abstract model
    V1 = VectorFunctionSpace('V1', domain, kind='H1')
    V2 = ScalarFunctionSpace('V2', domain, kind='L2')
    X  = ProductSpace(V1, V2)

    u0, u, v, du = elements_of(V1, names='u0, u, v, du')
    p0, p, q, dp = elements_of(V2, names='p0, p, q, dp')

    x, y, z  = domain.coordinates
    int_0 = lambda expr: integral(domain , expr)
    int_1 = lambda expr: integral(domain , expr)

    # time step
    dt = Constant(name='dt')

    # Boundaries
    bc = [EssentialBC(du, 0, boundary_h)]

    # Reynolds number
    Re = 0.004
    f  = Tuple(3,0,0)
    Fl = lambda u,p: Re**-1*inner(grad(u), grad(v)) - div(u)*q - p*div(v) + 1e-10*p*q
    F  = lambda u,p: dot(Transpose(grad(u))*u,v) + Fl(u,p)

    l = LinearForm((v, q), int_0(dot(u,v)-dot(u0,v) + dt/2 * (F(u,p) + F(u0,p0)) ) -int_1(dt*dot(f,v)))
    a = linearize(l, (u,p), trials=(du, dp))

    equation  = find((du, dp), forall=(v, q), lhs=a((du, dp), (v, q)), rhs=l(v, q), bc=bc)

    # Define (abstract) norms
    l2norm_du  = Norm(Matrix([du[0],du[1],du[2]]), domain, kind='l2')
    l2norm_dp  = Norm(dp     , domain, kind='l2')

    # ... create the computational domain from a topological domain
    domain_h = discretize(domain, ncells=ncells, comm=comm)

    # ... discrete spaces
    V1h = discretize(V1, domain_h, degree=degree, multiplicity=multiplicity)
    V2h = discretize(V2, domain_h, degree=degree, multiplicity=multiplicity)
    Xh  = V1h*V2h

    # ... discretize the equations
    equation_h = discretize(equation, domain_h, [Xh, Xh], backend=PSYDAC_BACKEND_GPYCCEL)

    a_h = equation_h.lhs
    l_h = equation_h.rhs

    # Discretize the norms
    l2norm_du_h = discretize(l2norm_du, domain_h, V1h, backend=PSYDAC_BACKEND_GPYCCEL)
    l2norm_dp_h = discretize(l2norm_dp, domain_h, V2h, backend=PSYDAC_BACKEND_GPYCCEL)

    u0_h = FemField(V1h)
    p0_h = FemField(V2h)

    u_h  = FemField(V1h)
    p_h  = FemField(V2h)

    du_h = FemField(V1h)
    dp_h = FemField(V2h)

    if store:
        output_m = OutputManager('nvs_spaces.yml', 'nvs_fields.h5')
        output_m.add_spaces(V1=V1h, V2=V2h)
        output_m.export_space_info()

    Tf = dt_h*(nt+1)
    t  = 0

    solver = petsc_solver if petsc else psydac_solver

    while t<Tf:
        t += dt_h
        print()
        print('======= time {}/{} ======='.format(t,Tf))

        u0_h[0].coeffs[:,:] = u_h[0].coeffs[:,:,:]
        u0_h[1].coeffs[:,:] = u_h[1].coeffs[:,:,:]
        u0_h[2].coeffs[:,:] = u_h[2].coeffs[:,:,:]
        p0_h.coeffs[:,:]    = p_h.coeffs[:,:,:]

        # Newton iteration
        for n in range(max_newton_iter):
            print()
            print('==== iteration {} ===='.format(n))

            M = a_h.assemble(u=u_h, p=p_h, dt=dt_h)
            b = l_h.assemble(u=u_h, p=p_h, u0=u0_h, p0=p0_h, dt=dt_h)

            apply_essential_bc(M, *equation_h.bc, identity=True)
            apply_essential_bc(b, *equation_h.bc)

            x,info = solver(M, b)

            du_h[0].coeffs[:] = x[0][:]
            du_h[1].coeffs[:] = x[1][:]
            du_h[2].coeffs[:] = x[2][:]
            dp_h.coeffs[:]    = x[3][:]

            # Compute L2 norm of increment
            l2_error_du = l2norm_du_h.assemble(du=du_h)
            l2_error_dp = l2norm_dp_h.assemble(dp=dp_h)

            print('L2_error_norm(du) = {}'.format(l2_error_du))
            print('L2_error_norm(dp) = {}'.format(l2_error_dp))

            if abs(l2_error_du) <= newton_tol:
                print()
                print('CONVERGED')
                break
            elif n == max_newton_iter-1 or abs(l2_error_du)>1/newton_tol or abs(l2_error_dp) > 1/newton_tol:
                print()
                print('NOT CONVERGED')
                t = Tf
                return solutions, p_h, domain, domain_h

            # update field
            u_h -= du_h
            p_h -= dp_h

        if store:
            output_m.add_snapshot(t=t, ts=t//dt_h-1)
            output_m.export_fields(u=u_h, p=p_h)

    if store:output_m.close()

def parse_input_arguments():

    import argparse

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description     = "Solve the Time dependent Navier Stokes equation on a 3D domain with" +
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

    parser.add_argument( '-ns',
        type    = int,
        default = 1,
        dest    = 'NSTEPS',
        metavar = 'NSTEPS',
        help    = 'Number of time-steps to be taken'
    )

    parser.add_argument( '-T',
        type    = float,
        dest    = 'END_TIME',
        metavar = 'END_TIME',
        help    = 'Run simulation until given final time'
    )
    # ...

    parser.add_argument( '-s',
        action  = 'store_true',
        dest    = 'store',
        help    = 'Save output files'
    )

    return parser.parse_args()
#------------------------------------------------------------------------------
if __name__ == '__main__':

    from time       import time
    from psydac.api.postprocessing import PostProcessManager

    Tf       = 1
    dt_h     = 0.01
    nt       = Tf//dt_h

    # Define topological domain 
    r_in  = 2.0 
    r_out = 4.0 
    A       = Cube('A', bounds1=(0, 18), bounds2=(r_in, r_out), bounds3=(0, 2*np.pi)) 
    mapp    = ConcentricPipe('M', 3, c1=0, c2=0, rmin=0, rmax=1) 
    domain  = mapp(A)

    boundary_h = Union(*[domain.get_boundary(axis=1, ext=-1), domain.get_boundary(axis=1, ext=1)])
    boundary_i = domain.get_boundary(axis=0, ext=-1)

    args = vars(parse_input_arguments())

    ncells       = args['ncells']
    degree       = args['degree']
    Tf           = args['END_TIME']
    nt           = args['NSTEPS']
    dt_h         = Tf/nt
    multiplicity = [1,1,1]

    store = args.get('store', None)

    run_time_dependent_navier_stokes_2d(domain, ncells, degree, multiplicity, boundary_h, boundary_i, dt_h=dt_h, nt=nt, petsc=False, store=store)

    if store:
        Pm = PostProcessManager(
            domain=domain,
            space_file='nvs_spaces.yml',
            fields_file='nvs_fields.h5',
            comm=comm,    
        )

        Pm.export_to_vtk(
            'visu_nvs',
            npts_per_cell=3,
            snapshots='all',
            fields=('u','p')
        )

        Pm.close()
