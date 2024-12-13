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
from MagInt.HubbardI_solver import Solver as Solver
from MagInt.Read_input import read_input_file
from MagInt.utils import *
import triqs.utility.mpi as mpi
from h5 import *
import numpy as np
import sys
import os


def run_shift_mu(mu_shift):
    """
    Applies a shift to of the chemical potential.

    This function adjusts the chemical potential by a specified amount, recalculates the atomic Self Energy and Green's
    function and stores the updated results in an HDF5 file.

    Parameters:
    -----------
    mu_shift : float
        The value by which to shift the chemical potential.

    Raises:
    -------
    FileNotFoundError
        If the input file 'maginteract.ini' is not found.
    RuntimeError
        If previous iterations are required but not found in the HDF5 file.
    Exception
        If there is an error reading previous data from the HDF5 file.

    Notes:
    ------
    - Reads parameters from 'maginteract.ini' and broadcasts them across MPI processes.
    - Initializes the `SumkDFT` object and solvers for each inequivalent shell.
    - Adjusts the chemical potential by `mu_shift` and recalculates Green's functions, self-energies,
      and effective atomic levels.
    - Stores updated results, including impurity self-energy and Green's functions, in the HDF5 file.

    """

    input_filename = 'maginteract.ini'
    if not os.path.isfile(input_filename):
        raise FileNotFoundError(f'Could not find input file {input_filename}.')

    print_logo()
    mpi.report('Reading the input file ' + input_filename)
    general_par, solver_par, basis_par, magint_par = read_input_file(input_filename)
    general_par['filename'] = general_par['folder'] + general_par['dft_filename']

    general_par = mpi.bcast(general_par)
    solver_par = mpi.bcast(solver_par)

    mpi.report('-' * 40)
    mpi.report('Starting Shift Chemical Potential calculation ... ')
    mpi.report('-' * 40)

    # Init the SumK class
    SK = SumkDFT(hdf_file=general_par['filename'] + '.h5', use_dft_blocks=False)

    # Initialization U
    U_full = []
    U_dc = []
    J_dc = []
    useval = []

    nlm = []
    nlms = []
    mpi.report('-' * 40)
    for i_sh in range(SK.n_inequiv_shells):
        nlm.append(2 * solver_par['ls'][i_sh] + 1)
        nlms.append(2 * nlm[i_sh])
        # Initialization U
        if solver_par['kanamori'] is True:
            U_full_val, U_dc_val, J_dc_val = Kanamori_hamiltonian(U_int=solver_par['u_int'][i_sh],
                                                                  J_hund=solver_par['j_hund'][i_sh],
                                                                  orbital=solver_par['ls'][i_sh],
                                                                  dft_exec=general_par['dft_exec'])
            U_full.append(U_full_val)
            if mpi.is_master_node():
                print('\nKanamori U from U_int : ', solver_par['u_int'][i_sh], ' and J_hund : ',
                      solver_par['j_hund'][i_sh],
                      'for shell ', i_sh)
        else:
            U_dc_val = solver_par['u_int'][i_sh]
            J_dc_val = solver_par['j_hund'][i_sh]
            if mpi.is_master_node():
                print('\nFull U from U_int : ', solver_par['u_int'][i_sh], ' and J_hund : ', solver_par['j_hund'][i_sh],
                      'for shell ', i_sh)

        U_dc.append(U_dc_val)
        J_dc.append(J_dc_val)

        # Initialization DC
        useval.append(
            U_dc[i_sh] * (solver_par['n_el'][i_sh] - 0.5) - J_dc[i_sh] * (solver_par['n_el'][i_sh] * 0.5 - 0.5))

    # check if there are previous runs:
    previous_present = False
    if mpi.is_master_node():
        f = HDFArchive(general_par['filename'] + '.h5', 'a')
        if 'dmft_output' in f:
            print('\nPrevious iterations are present!')
            ar = f['dmft_output']
            if 'iterations' in ar:
                previous_present = True
        else:
            del f
            raise RuntimeError(
                "Previous iterations are not found! Please ensure that in the .h5 file a 'dmft_output' "
                "group is present")
        del f
    previous_present = mpi.bcast(previous_present)

    S = []
    dm = []
    for i_sh in range(SK.n_inequiv_shells):
        S.append(Solver(beta=solver_par['beta'], l=solver_par['ls'][i_sh], n_iomega=solver_par['n_iomega'],
                        Nmoments=solver_par['n_moments'], use_spin_orbit=solver_par['use_spin_orbit']))

        if (previous_present):
            try:
                # load previous data:
                mpi.report("Using stored data for initialisation")
                if mpi.is_master_node():
                    ar = HDFArchive(general_par['filename'] + '.h5', 'a')
                    S[i_sh].Sigma_iw << ar['dmft_output']['Sigma_iw_' + str(i_sh)]
                    del ar
                    SK.chemical_potential, SK.dc_imp, SK.dc_energ = SK.load(
                        ['chemical_potential', 'dc_imp', 'dc_energ'])
                    print('\nChemical  potential : ', SK.chemical_potential, '\n')
                S[i_sh].Sigma_iw << mpi.bcast(S[i_sh].Sigma_iw)
                SK.chemical_potential = mpi.bcast(SK.chemical_potential)
                SK.dc_imp = mpi.bcast(SK.dc_imp)
                SK.dc_energ = mpi.bcast(SK.dc_energ)
            except Exception as e:
                raise RuntimeError(f"Error in reading old data from h5: {e}")

        # calculated fixed DC at the beginning
        dm.append(S[i_sh].G_iw.density())
        SK.calc_dc(dm[i_sh], U_interact=solver_par['u_int'][i_sh], J_hund=solver_par['j_hund'][i_sh], orb=i_sh,
                   use_dc_formula=solver_par['dc_type'], use_dc_value=useval[i_sh])
        if not previous_present:
            S[i_sh].Sigma_iw << SK.dc_imp[i_sh]['ud'][0, 0]

    # Calculate Sigma with shifted mu
    if mpi.is_master_node():
        print('-' * 40)
        print('\nOld chemical potential : ', SK.chemical_potential, '\n')
    SK.chemical_potential = SK.chemical_potential + mu_shift
    if mpi.is_master_node():
        print('New chemical potential : ', SK.chemical_potential, '\n')

    # put Sigma into the SumK class:
    SK.set_Sigma([S[i].Sigma_iw for i in range(SK.n_inequiv_shells)])

    for i_sh in range(SK.n_inequiv_shells):
        mpi.report('\n\nStarting calculation for shell ' + str(i_sh) + '\n\n')

        S[i_sh].G_iw <<= SK.extract_G_loc()[i_sh]
        mpi.report("Total charge of Gloc : ", S[i_sh].G_iw.total_density().real)

        eal = SK.eff_atomic_levels()
        for s in eal[i_sh]:
            eal[i_sh][s].real[abs(eal[i_sh][s].real) < solver_par['tol_eal']] = 0
            eal[i_sh][s].imag[abs(eal[i_sh][s].imag) < solver_par['tol_eal']] = 0

        S[i_sh].set_atomic_levels(eal=eal[i_sh])

        if solver_par['kanamori'] is False:
            S[i_sh].solve(U_int=solver_par['u_int'][i_sh], J_hund=solver_par['j_hund'][i_sh],
                          verbosity=general_par['verbosity'])
        else:
            S[i_sh].solve(u4ind=U_full[i_sh], verbosity=general_par['verbosity'])

        mpi.barrier()
        if mpi.is_master_node():
            os.rename('eal.dat', 'eal_site_%s.dat' % i_sh)
            os.rename('atomic_levels.dat', 'atomic_levels_site_%s.dat' % i_sh)
            os.rename('states.dat', 'states_site_%s.dat' % i_sh)

        # store the impurity self-energy, GF as well as correlation energy in h5
        if mpi.is_master_node():
            ar = HDFArchive(general_par['filename'] + '.h5', 'a')
            try:
                ar['dmft_output']['G_iw_' + str(i_sh)] = S[i_sh].G_iw
                ar['dmft_output']['Sigma_iw_' + str(i_sh)] = S[i_sh].Sigma_iw
            except:
                pass
            del ar

        SK.save(['chemical_potential', 'dc_imp', 'dc_energ'])
