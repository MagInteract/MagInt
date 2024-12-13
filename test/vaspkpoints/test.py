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
from MagInt.utils import compare_lists_arrays
from MagInt.Read_kp_sym import *
from h5 import *
import numpy as np
import os

dft_filename = os.getcwd().rpartition('/')[2]

general_par = {}
general_par['verbosity'] = 2
general_par['filename'] = 'vasp'
general_par['folder'] = '.'
general_par['dft_exec'] = 'Vasp'
general_par['use_symmetries'] = True

KP = read_kp_sym(general_par)

h = HDFArchive('ref.h5', 'r')

assert np.allclose(KP.nk, h['KP']['KP_nk'])
assert np.allclose(KP.weight, h['KP']['KP_weight'])


