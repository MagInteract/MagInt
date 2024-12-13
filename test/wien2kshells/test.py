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

from MagInt.utils import compare_lists_arrays, compare_dicts, compare_lists
from MagInt.Shells import *
import numpy
from h5 import *


general_par = {}
general_par['dft_exec'] = 'Wien2k'
general_par['filename'] = 'case'
general_par['verbosity'] = 2
general_par['tol_shells'] = 0.0001

Interactions = (('Re10', 'Re10'),)
interact_sites = {}
interact_sites['Re10'] = [0]
interact_basis = 1 

int_pairs = {}

for int0 in Interactions:
    if general_par['verbosity'] > 0:
        mpi.report("\n----INTERACTION TYPE: %s-%s" % (int0[0], int0[1]))
    int_pairs[int0] = shells(general_par, 'case', int0, interact_sites, interact_basis, max_n_shells=2, print_vecs='Basis')

#if (mpi.is_master_node()):
#    ar = HDFArchive('ref.h5', 'a')
#    if not ('struct' in ar):
#        ar.create_group('struct')
#    # data on types, interacting pairs, ground-state multiplicity
#    ar['struct']['el_names'] = int_pairs[int0].struct.el_names
#    ar['struct']['q_types'] = int_pairs[int0].struct.q_types
#    ar['struct']['nions'] = int_pairs[int0].struct.nions
#    ar['struct']['type_of_ion'] = int_pairs[int0].struct.type_of_ion
#    ar['struct']['a_brav'] = int_pairs[int0].struct.a_brav
   

h = HDFArchive('ref.h5', 'r')

assert compare_lists_arrays(int_pairs[int0].struct.q_types, h['struct']['q_types'])
assert compare_lists(int_pairs[int0].struct.type_of_ion, h['struct']['type_of_ion'])
assert compare_lists(int_pairs[int0].struct.nions, h['struct']['nions'])
assert np.allclose(int_pairs[int0].struct.a_brav, h['struct']['a_brav'])


