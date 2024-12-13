# ------------------------------------------------------------------------------------#
# MagInteract library
# Written by  : Leonid V. Pourovskii (CPHT Ecole Polytechnique) 2012-2023
#             : Dario Fiore Mosca (Ecole Polytechnique) 2023
# Email: leonid@cpht.polytechnique.fr
# ------------------------------------------------------------------------------------#
#
#    MagInteract implements the force-theorem Hubbard-I (FT-HI) approach to intersite
#    exchange interactions in correlated insulators.
#    The formalism is given in L. V. Pourovskii Phys. Rev. B 94, 115117 (2016)
#    This Python-3 version is based on TRIQS library
#
# ------------------------------------------------------------------------------------#

from h5 import *
from MagInt.HubbardI_solver import Solver
from triqs.utility.comparison_tests import *
import numpy

l = 2
nlm = 2 * l + 1

S = Solver(beta=200, l=l)

eal = {}
eal['up'] = -1 * numpy.identity(5)
eal['down'] = -1 * numpy.identity(5)
S.set_atomic_levels(eal=eal)

S.solve(U_int=6.0, J_hund=0.6)

h = HDFArchive('ref.h5', 'r')
assert_block_gfs_are_close(S.G_iw, h['G'])
assert_block_gfs_are_close(S.Sigma_iw, h['Sigma'])
