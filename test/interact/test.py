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
from MagInt.HubbardI_interact import HubbardI_interact
from triqs.utility.comparison_tests import *
import triqs.utility.mpi as mpi
import numpy as np

l = 2

S = HubbardI_interact(beta=200, l=l, U_int=6.0, J_hund=0.6, n_lev=5, use_spin_orbit=True)

eal = {}
eal['ud'] = np.zeros((10, 10))
eal['ud'][:5, :5] = 1 * np.identity(5)
eal['ud'][5:10, 5:10] = -1 * np.identity(5)

S.set_ud_levels(eal=eal)
S.run_HI()

h = HDFArchive('ref.h5', 'r')
assert_block_gfs_are_close(S.G_iw, h['G_iw'])
assert_block_gfs_are_close(S.Sigma_iw, h['Sigma_iw'])
