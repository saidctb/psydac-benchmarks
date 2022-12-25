# -*- coding: UTF-8 -*-

import os
import numpy as np

from time import time
setup_time1 = time()
from sympde.topology import Mapping,IdentityMapping
from psydac.api.postprocessing import PostProcessManager
from psydac.api.postprocessing import OutputManager
from sympy import sin, pi

from sympde.calculus import laplace, dot, grad,div
from sympde.topology import Cube, Mapping, Domain, Union
from sympde.topology import ScalarFunctionSpace, elements_of
from sympde.topology import NormalVector
from sympde.expr     import BilinearForm, LinearForm, Norm, TerminalExpr
from sympde.expr     import EssentialBC, find
from sympde.expr     import integral

from psydac.api.discretization import discretize
from psydac.api.settings       import PSYDAC_BACKEND_GPYCCEL, PSYDAC_DEFAULT_FOLDER
from sympy import lambdify
from mpi4py import MPI

def remove_folder(path):
    os.system('rm -rf "%s" &' % path)

class SphericalMapping(Mapping):
    """
    Represents a Spherical 3D Mapping object.
    """
    _expressions = {'x': 'x1*sin(x2)*cos(x3)',
                    'y': 'x1*sin(x2)*sin(x3)',
                    'z': 'x1*cos(x2)'}

    _ldim        = 3
    _pdim        = 3

class PolarMapping(Mapping):
    """
    Represents a Polar 3D Mapping object (Annulus).

    Examples

    """
    _expressions = {'x': 'c1 + (rmin*(1-x1)+rmax*x1)*cos(x2)',
                    'y': 'c2 + (rmin*(1-x1)+rmax*x1)*sin(x2)',
                    'z': 'x3'}

    _ldim        = 3
    _pdim        = 3

def construct_mapping(ncells, degree, comm):

    from sympde.topology                    import Cube
    from psydac.cad.geometry                import Geometry
    from psydac.fem.splines                 import SplineSpace
    from psydac.fem.tensor                  import TensorFemSpace
    from psydac.mapping.discrete            import SplineMapping
    from psydac.ddm.cart                    import DomainDecomposition

    p1 , p2 , p3  = degree
    nc1, nc2, nc3 = ncells

    filename = 'quarter_annulus_3d_{}_{}.h5'.format(ncells[0], degree[0])
    rank = 0 if comm is None else comm.rank

    if rank == 0:
        # Create the domain decomposition
        domain_decomposition = DomainDecomposition(ncells=[nc1,nc2,nc3], periods=[False]*3)

        default_params = dict( rmin=0.0, rmax=1.0, c1=0.0, c2=0.0)
        r_in    = 1.0
        r_out   = 2.0
        lims1   = (r_in, r_out)
        lims2   = (np.pi/4, 3*np.pi/4)
        lims3   = (0, np.pi/2)
        domain = Cube('Omega', bounds1=lims1, bounds2=lims2, bounds3=lims3)

        # Create tensor spline space, distributed
        V1    = SplineSpace( grid=np.linspace( *lims1, num=nc1+1 ), degree=p1, periodic=False )
        V2    = SplineSpace( grid=np.linspace( *lims2, num=nc2+1 ), degree=p2, periodic=False )
        V3    = SplineSpace( grid=np.linspace( *lims3, num=nc3+1 ), degree=p3, periodic=False )
        space = TensorFemSpace( domain_decomposition, V1, V2, V3 )

        map_analytic = SphericalMapping( 'M', dim=3 )
        # Create spline mapping by interpolating analytical one

        # Topological domain
        domain  = map_analytic(domain)
        mapping = SplineMapping.from_mapping( space, map_analytic )

        # Define ncells as a dict
        mappings = {domain.name: mapping}
        ncells   = {domain.name:[len(space.breaks)-1 for space in mapping.space.spaces]}
        periodic = {domain.name:[space.periodic for space in mapping.space.spaces]}

        # create a geometry from a topological domain and the dict of mappings
        geo = Geometry(domain=domain, ncells=ncells, periodic=periodic, mappings=mappings)

        # Export to file
        geo.export(filename)

#==============================================================================
def run_model(filename, ncells, degree, comm, backend):

    backend['folder'] = "poisson_3d_psydac_{}_{}_{}_{}_{}".format(ncells[0], degree[0], comm.size, int(os.environ.get('OMP_NUM_THREADS', 1)), filename is None)
    backend['flags']  = "-O3 -march=native -mtune=native  -mavx -ffast-math -ffree-line-length-none"
    PSYDAC_DEFAULT_FOLDER['name'] = '__psydac__' + backend['folder']

    # Define topological domain
    Omega = Domain.from_file(filename=filename)

    # Method of manufactured solutions: define exact
    # solution u_e, then compute right-hand side f
    x, y, z = Omega.coordinates

    u_e   = sin(pi*x)*sin(pi*y)*sin(pi*z)
    f     = laplace(laplace(u_e))
    f     = TerminalExpr(f, Omega)

    # Define abstract model
    V = ScalarFunctionSpace('V', Omega)
    u, v = elements_of(V, names='u, v')

    nn = NormalVector('nn')

    kappa = 10*ncells[0]*degree[0]
    expr_b = - laplace(u)*dot(grad(v), nn)\
             - dot(grad(u), nn)*laplace(v) \
             + kappa*dot(grad(u),nn)*dot(grad(v),nn)

    a = BilinearForm((u,v), integral(Omega, laplace(v) * laplace(u)) + integral(Omega.boundary, expr_b))

    expr_b = - dot(grad(u_e), nn)*laplace(v) \
             + kappa*dot(grad(u_e),nn)*dot(grad(v),nn)

    l =   LinearForm(   v , integral(Omega, f * v) + integral(Omega.boundary, expr_b))

    bc = [EssentialBC(  u, u_e, Omega.boundary)]
#         EssentialBC(dot(grad(u), nn), dot(grad(u_e), nn), Omega.boundary)]

    equation = find(u, forall=v, lhs=a(u,v), rhs=l(v), bc=bc)

    # Define (abstract) error norms
    error  = u - u_e
    l2norm = Norm(error, Omega, kind='l2')
    h1norm = Norm(error, Omega, kind='h1')
    h2norm = Norm(error, Omega, kind='h2')

    # Create computational domain from topological domain
    Omega_h = discretize(Omega, filename=filename, comm=comm)

    # Create discrete spline space
    Vh = discretize(V, Omega_h,degree=degree)

    # Discretize equation
    equation_h = discretize(equation, Omega_h, [Vh, Vh], backend=PSYDAC_BACKEND_GPYCCEL)

    # Discretize norms
    l2norm_h = discretize(l2norm, Omega_h, Vh, backend=PSYDAC_BACKEND_GPYCCEL)
    h1norm_h = discretize(h1norm, Omega_h, Vh, backend=PSYDAC_BACKEND_GPYCCEL)
    h2norm_h = discretize(h2norm, Omega_h, Vh, backend=PSYDAC_BACKEND_GPYCCEL)

    # Solve discrete equation to obtain finite element coefficients
    equation_h.set_solver(solver='cg' ,tol=1e-18, maxiter=100, info=True, verbose=False)

    comm.Barrier()
    try:
        remove_folder(backend['folder'])
        remove_folder(PSYDAC_DEFAULT_FOLDER['name'])
    except:
        pass
    comm.Barrier()

    setup_time2 = time()
    T = comm.reduce(setup_time2-setup_time1,op=MPI.MAX)

    infos = {}
    infos['title'] = 'biharmonic_3d'
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
    b  = lhs.assemble()
    t2 = time()
    T = comm.reduce(t2-t1,op=MPI.MAX)

    infos['blinear_form_assembly_time2'] = T
    b   = rhs.assemble()
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

    equation_h.set_solver('cg', tol=1e-18, maxiter=100, info=True)
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
        name = (infos['title'],) + (('geof',) if filename else ()) + infos['ncells'] + infos['degree'] + (comm.size, infos['number_of_threads'])
        name = '_'.join([str(i) for i in name])
        np.save('results/' + name, infos)
    
    u_h,info = equation_h.solve()


    return locals()

#==============================================================================
def parse_input_arguments():

    import argparse

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description     = "Solve the biharmonic equation on a 3D domain with" +
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


    parser.add_argument('-a', action='store_true', \
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

#==============================================================================
def main(degree, ncells, **kwargs):

    comm = MPI.COMM_WORLD
    rank = comm.rank

#    construct_mapping(ncells, degree, comm)
    filename = 'quarter_annulus_3d_{}_{}.h5'.format(ncells[0], degree[0])
    namespace = run_model(filename, ncells, degree, comm, backend=PSYDAC_BACKEND_GPYCCEL)

    return namespace
#==============================================================================
if __name__ == '__main__':

    args = parse_input_arguments()
    args = vars(args)
    namespace = main( **args )

