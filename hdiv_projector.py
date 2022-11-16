import numpy as np

from sympde.topology import Cube
from sympde.topology import VectorFunctionSpace

from psydac.feec.global_projectors import Projector_Hdiv
from psydac.api.settings           import PSYDAC_BACKEND_GPYCCEL
from psydac.api.discretization     import discretize
from psydac.linalg.block           import BlockVector

from time import time
from mpi4py import MPI

def run_hdiv_projector(ncells, degree, comm):

    backend['folder'] = "hdiv_psydac_{}_{}_{}".format(ncells[0], degree[0], comm.size)
    backend['flags']  = "-O3 -march=native -mtune=native  -mavx -ffast-math"

    #+++++++++++++++++++++++++++++++
    # 1. Abstract model
    #+++++++++++++++++++++++++++++++
    V  = VectorFunctionSpace('V', domain, kind='hdiv')

    domain_h = discretize(domain, ncells=ncells, comm=comm)
    Vh       = discretize(V, domain_h, degree=degree, basis='M')

    infos = {}
    infos['title'] = 'hdiv_projector_3d'
    infos['ncells'] = tuple(ncells)
    infos['degree'] = tuple(degree)
 
    comm.Barrier()
    t0 = time()
    P2 = Projector_Hdiv(Vh)
    b  = BlockVector(Vh.vector_space)
    t1 = time()

    tt = comm.reduce(t1-t0,op=MPI.MAX)

    infos['matrix_assembly'] = tt

    comm.Barrier()
    t0 = time()
    P2._solver.solve(b)
    t1 = time()

    tt = comm.reduce(t1-t0,op=MPI.MAX)

    infos['solver'] = tt

    if comm.rank == 0:
        name = (infos['title'],) + infos['ncells'] + infos['degree'] + (comm.size,)
        name = '_'.join([str(i) for i in name])
        np.save('results/' + name, infos)

def parse_input_arguments():

    import argparse

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description     = "Run Hdiv projector"
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

    return parser.parse_args()
if __name__ == '__main__':

    domain = Cube('domain',bounds1=(0, 1), bounds2=(0., 1), bounds3=(0., 1))

    x,y,z = domain.coordinates

    args = vars(parse_input_arguments())

    ncells = args['ncells']
    degree = args['degree']

    comm    = MPI.COMM_WORLD
    backend = PSYDAC_BACKEND_GPYCCEL

    run_hdiv_projector(ncells, degree, comm)
