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
from psydac.api.settings             import PSYDAC_BACKEND_GPYCCEL
from psydac.feec.pull_push           import pull_2d_hcurl
from psydac.fem.vector               import ProductFemSpace
from psydac.linalg.iterative_solvers import pcg
from psydac.api.postprocessing       import OutputManager
from psydac.linalg.utilities         import array_to_psydac

def remove_folder(path):
    shutil.rmtree(path)
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

def splitting_integrator_sp(e0, bhalf, M1, M2, R1, R2, S1, S2, R1_b, BC, J, dt, niter, Vh, l2norm_h):
    from scipy.sparse.linalg import splu
    e_history = [e0.toarray()]
    b_history = [bhalf.toarray()]

    M1 = splu((M1).tosparse().tocsc())
    M2 = splu(M2.tosparse().tocsc())
    R1 = R1.tosparse().tocsr()
    R2 = R1
    S1 = S1.tosparse().tocsr()
    S2 = S2.tosparse().tocsr()
    R1_b = R1_b.tosparse().tocsr()

    t = 0.
    for ts in range(niter):
        e = e_history[ts]
        b = b_history[ts]

        y1    = dt * R1.dot(b)  -dt*R1_b.dot(b) +dt*(J(t=t+dt/2).toarray()) + dt * S1.dot(b)
        e_new = e + M1.solve(y1)

        y2    = -dt * R2.dot(e_new) - dt * S2.dot(e_new) +dt*BC(t=t+dt).toarray()
        b_new = b + M2.solve(y2)
        e_history.append(e_new)
        b_history.append(b_new)

        t += dt

    for i,(e,b) in enumerate(zip(e_history,b_history)):
        e_history[i] = array_to_psydac(e, Vh.vector_space)
        b_history[i] = array_to_psydac(b, Vh.vector_space)

    return e_history, b_history

#===============================================================================
def splitting_integrator(e0, bhalf, M1, M2, R1, R2, dt, niter, Vh, l2norm_h,cg_niter):

    e_history = [e0]
    b_history = [bhalf]

    t = 0.
    for ts in range(niter):
        e = e_history[ts]
        b = b_history[ts]

        y1    = dt*(R1.dot(b))
        e_new = e + cg(M1, y1, tol=1e-15, maxiter=cg_niter)[0]

        y2    = dt * (R2.dot(e_new))
        b_new = b + cg(M2, y2, tol=1e-15, maxiter=cg_niter)[0]

        e_history.append(e_new)
        b_history.append(b_new)

        t += dt

    return e_history, b_history
#==============================================================================
def run_maxwell_3d(Eex, Bex, J, domain, ncells, degree,  dt, niter, T, backend, comm, cg_niter,store=None):

    backend['folder'] = "maxwell_3d_psydac_{}_{}_{}".format(ncells[0], degree[0], comm.size)
    backend['flags']  = "-O3 -march=native -mtune=native  -mavx -ffast-math"

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
    a2_h = discretize(a2, domain_h, [Vh, Vh], backend=backend)
    a3_h = discretize(a3, domain_h, [Vh, Vh], backend=backend)

    l1_h = discretize(l1, domain_h, Vh, backend=backend)
    l2_h = discretize(l2, domain_h, Vh, backend=backend)
    l3_h = discretize(l3, domain_h, Vh, backend=backend)
    l4_h = discretize(l4, domain_h, Vh, backend=backend)

    l2norm_h   = discretize(l2norm, domain_h, Vh, backend=backend)

    comm.Barrier()
    if comm.rank == 0:
        try:
            remove_folder(backend['folder'])
        except:
            pass
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
    M2 = M1
    R1 = a2_h.assemble()
    R2 = a3_h.assemble()
    t2 = time()

    tt = comm.reduce(t2-t1,op=MPI.MAX)

    # store bilinear forms assembly timing
    infos['bilinear_form_assembly_time'] = tt

    e0    = cg(M1, l1_h.assemble(t=0.), tol=1e-1, maxiter=10)[0]
    bhalf = cg(M2, l2_h.assemble(t=dt/2), tol=1e-1, maxiter=10)[0]

    b   = l1_h.assemble(t=0.)
    out = b.copy()
    st  = 0
    for i in range(20):
        comm.Barrier()
        t1 = time()
        M1.dot(b, out=out)
        t2 = time()
        b.ghost_regions_in_sync = False
        st += t2-t1

    tt = comm.reduce(st/20,op=MPI.MAX)

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

    comm.Barrier()
    t0       = time()
    e_history, b_history = splitting_integrator(e0, bhalf, M1, M2, R1, R2, dt, niter, Vh, l2norm_h, cg_niter)
    t1       = time()

    tt = comm.reduce(t1-t0,op=MPI.MAX)

    # store solution timing
    infos['solution'] = tt

    Ef = FemField(Vh, coeffs=e_history[-1])
    l2_error = l2norm_h.assemble(E=Ef, t=T)

    infos['l2_error'] = l2_error

    if comm.rank == 0:
        name = (infos['title'],) + infos['ncells'] + infos['degree'] + (comm.size, infos['number_of_threads'])
        name = '_'.join([str(i) for i in name])
        np.save('results/' + name, infos)

    if store:
        output_m = OutputManager('maxwell_space.yml', 'maxwell_fields.h5')
        output_m.add_spaces(V=Vh)
        output_m.export_space_info() # Writes the space information to Yaml

        for i,e in enumerate(e_history):
            e = FemField(Vh, coeffs=e)
            output_m.add_snapshot(t=i*dt, ts=i)
            output_m.export_fields(e=e)

        output_m.close()

    return infos, l2_error

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

    if store:
        Pm = PostProcessManager(
            domain=domain,
            space_file='maxwell_space.yml',
            fields_file='maxwell_fields.h5',
        )

        Pm.export_to_vtk(
            'visu_maxwell',
            npts_per_cell=3,
            snapshots='all',
            fields=('e')
        )

        Pm.close()

