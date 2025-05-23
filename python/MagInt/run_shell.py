# ------------------------------------------------------------------------------------#
# MagInt library
# Written by  : Leonid V. Pourovskii (CPHT Ecole Polytechnique) 2012-2024
#             : Dario Fiore Mosca (CPHT Ecole Polytechnique) 2023-2024
# Email: leonid@cpht.polytechnique.fr
# ------------------------------------------------------------------------------------#
#
#    MagInt implements the force-theorem Hubbard-I (FT-HI) approach to intersite
#    exchange interactions in correlated insulators.
#    The formalism is given in L. V. Pourovskii Phys. Rev. B 94, 115117 (2016)
#    This Python-3 version is based on TRIQS library
#
# ------------------------------------------------------------------------------------#
from triqs_dft_tools.sumk_dft import *
from triqs_dft_tools.converters.wien2k import *
from MagInt.Shells import *
from MagInt.Read_input import read_input_file
from MagInt.utils import *
import triqs.utility.mpi as mpi
from h5 import *
import numpy as np
import sys
import os


def run_shell():
    """
    Perform shell-based calculations for magnetic interactions.

    This function reads input parameters, initializes interaction types, sites, and basis, and performs shell evaluations
    for specified interactions.

    Parameters:
    -----------
    None

    Raises:
    -------
    FileNotFoundError
        If the input file 'maginteract.ini' is not found.

    Notes:
    ------
    - Reads parameters from 'maginteract.ini' and broadcasts them across MPI processes.
    - Handles initialization of interaction types, sites, and basis.
    - If the DFT execution is Wien2k, lattice parameters are read from the output file '.outputnn'.
    - Performs shell evaluation for specified interactions using the `shells` function.

    """

    input_filename = 'maginteract.ini'
    if not os.path.isfile(input_filename):
        raise FileNotFoundError(f'Could not find input file {input_filename}.')

    mpi.report('Reading the input file ' + input_filename)
    general_par, solver_par, basis_par, magint_par = read_input_file(input_filename)
    general_par['filename'] = general_par['folder'] + general_par['dft_filename']

    fname = general_par['filename']
    general_par = mpi.bcast(general_par)

    interact_types = {}  # set()
    interact_sites = {}
    interact_basis = {}

    if general_par['dft_exec'] == 'Wien2k':
        # latt = np.zeros([3])
        f = open(fname + ".outputnn")
        head = [f.readline() for i in range(4)]
        width = 10
        #rela_line = head[3].strip()
        rela_line=head[3]
        latt = np.array([float(rela_line[j:j + width].strip()) for j in range(0, 3 * width, width)])
        if general_par['verbosity'] > 0:
            mpi.report("\nA = %10.6f    B = %10.6f    C = %10.6f" % (latt[0], latt[1], latt[2]))
        f.close()

        Interactions = (('Ion0', 'Ion1'),)

        for int in Interactions:
            for type in int:
                # self.interact_types.add(type)
                interact_types[type] = None
                if type not in interact_sites:
                    interact_sites[type] = []
                    interact_basis[type] = 0
                    count = 0
                    for icrsh in range(1):
                        interact_sites[type].append((icrsh, count))
                        interact_basis[type] += 1
                        count += 1

        for int in Interactions:
            mpi.report('-' * 40)
            mpi.report('Starting Shell evaluation ... ')
            mpi.report('-' * 40)
            shells(general_par, fname, int, interact_sites, interact_basis, max_n_shells=3, run_only=True)
