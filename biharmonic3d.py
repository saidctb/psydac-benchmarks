# -*- coding: UTF-8 -*-

import numpy as np

from sympde.topology import Mapping,IdentityMapping
from psydac.api.postprocessing import PostProcessManager
from psydac.api.postprocessing import OutputManager

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

def construct_mapping(ncells, degree, comm=None):

    from sympde.topology                    import Cube
    from psydac.cad.geometry                import Geometry
    from psydac.fem.splines                 import SplineSpace
    from psydac.fem.tensor                  import TensorFemSpace
    from psydac.mapping.discrete            import SplineMapping
    from psydac.ddm.cart                    import DomainDecomposition

    p1 , p2 , p3  = degree
    nc1, nc2, nc3 = ncells

    rank = 0 if comm is None else comm.rank

    if rank == 0:
        # Create the domain decomposition
        domain_decomposition = DomainDecomposition(ncells=[nc1,nc2,nc3], periods=[False]*3)

        default_params = dict( rmin=0.0, rmax=1.0, c1=0.0, c2=0.0)
        r_in    = 1.0
        r_out   = 4.0
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
        geo.export('quarter_annulus_3d.h5')

#==============================================================================
def run_model(ncells, degree, comm=None, store=None):

    from sympy import sin, pi

    from sympde.calculus import laplace, dot, grad,div
    from sympde.topology import Cube, Mapping, Domain, Union
    from sympde.topology import ScalarFunctionSpace, elements_of
    from sympde.topology import NormalVector
    from sympde.expr     import BilinearForm, LinearForm, Norm, TerminalExpr
    from sympde.expr     import EssentialBC, find
    from sympde.expr     import integral

    from psydac.api.discretization import discretize
    from psydac.api.settings       import PSYDAC_BACKEND_GPYCCEL

    # Define topological domain
    Omega = Domain.from_file(filename='quarter_annulus_3d.h5')

    # Method of manufactured solutions: define exact
    # solution u_e, then compute right-hand side f
    x, y, z = Omega.coordinates

    u_e   = (x**2+y**2+z**2)**2
    f     = laplace(laplace(u_e))
    f     = TerminalExpr(f, Omega)

    # Define abstract model
    V = ScalarFunctionSpace('V', Omega)
    u, v = elements_of(V, names='u, v')

    nn = NormalVector('nn')

    kappa = 10**3
    expr_b = - laplace(u)*dot(grad(v), nn)\
             - dot(grad(u), nn)*laplace(v) \
             + kappa*dot(grad(u),nn)*dot(grad(v),nn)

    a = BilinearForm((u,v), integral(Omega, laplace(v) * laplace(u)) + integral(Omega.boundary, expr_b))

    expr_b = - dot(grad(u_e), nn)*laplace(v) \
             + kappa*dot(grad(u_e),nn)*dot(grad(v),nn)

    l =   LinearForm(   v , integral(Omega, f * v) + integral(Omega.boundary, expr_b))

    bc = [EssentialBC(               u, u_e, Omega.boundary)]
#          EssentialBC(dot(grad(u), nn), dot(grad(u_e), nn), Omega.boundary)]

    equation = find(u, forall=v, lhs=a(u,v), rhs=l(v), bc=bc)

    # Define (abstract) error norms
    error  = u - u_e
    l2norm = Norm(error, Omega, kind='l2')
    h1norm = Norm(error, Omega, kind='h1')
    h2norm = Norm(error, Omega, kind='h2')

    # Create computational domain from topological domain
    Omega_h = discretize(Omega, filename='quarter_annulus_3d.h5', comm=comm)

    # Create discrete spline space
    Vh = discretize(V, Omega_h,degree=degree)

    # Discretize equation
    equation_h = discretize(equation, Omega_h, [Vh, Vh], backend=PSYDAC_BACKEND_GPYCCEL)

    # Discretize norms
    l2norm_h = discretize(l2norm, Omega_h, Vh, backend=PSYDAC_BACKEND_GPYCCEL)
    h1norm_h = discretize(h1norm, Omega_h, Vh, backend=PSYDAC_BACKEND_GPYCCEL)
    h2norm_h = discretize(h2norm, Omega_h, Vh, backend=PSYDAC_BACKEND_GPYCCEL)

    # Solve discrete equation to obtain finite element coefficients
    equation_h.set_solver(solver='cg' ,tol=1e-14, maxiter=10**10, info=True)
    u_h,info = equation_h.solve()

    # Compute error norms from solution field
    l2_error = l2norm_h.assemble(u=u_h)
    h1_error = h1norm_h.assemble(u=u_h)
    h2_error = h2norm_h.assemble(u=u_h)

    if store:
        output_m = OutputManager('biharmonic_spaces.yml', 'biharmonic_fields.h5')

        output_m.add_spaces(V=Vh)
        output_m.export_space_info()

        output_m.set_static()
        output_m.export_fields(u=u_h)
        output_m.close()
    return l2_error, h1_error, h2_error, info

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

    parser.add_argument( '-s',
        action  = 'store_true',
        dest    = 'store',
        help    = 'Save output files'
    )

    return parser.parse_args()

#==============================================================================
def main(degree, ncells, store=False):

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.rank

    construct_mapping(ncells, degree, comm)
    l2_error, h1_error, h2_error, info = run_model(ncells, degree, comm, store=store)

    if rank == 0:
        print()
        print('L2 error = {}'.format(l2_error))
        print('H1 error = {}'.format(h1_error))
        print('H2 error = {}'.format(h2_error))
        print('CG info  = {}'.format(info))
        print(flush=True)

    if comm:
        comm.Barrier()

#==============================================================================
if __name__ == '__main__':

    args = parse_input_arguments()
    args = vars(args)
    main( **args )

    if args.get('store',None):
        Pm = PostProcessManager(
            domain=domain,
            space_file='biharmonic_spaces.yml',
            fields_file='biharmonic_fields.h5',
            comm=comm,    
        )

        Pm.export_to_vtk(
            'visu_biharmonic',
            npts_per_cell=3,
            snapshots='all',
            fields=('u',)
        )

        Pm.close()

