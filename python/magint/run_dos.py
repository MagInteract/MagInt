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
from triqs_dft_tools.sumk_dft_tools import *
from triqs_dft_tools.converters.wien2k import *
from MagInt.HubbardI_solver import Solver as Solver
from MagInt.Read_input import read_input_file
import triqs.utility.mpi as mpi
from h5 import *
from MagInt.utils import *
import sys


def run_dos():
    """
    Perform a Density of States (DOS) calculation.

    This function initializes parameters, sets up solvers, handles previous iterations if available,
    and calculates the DOS using input from a configuration file and interaction parameters.

    Parameters:
    -----------
    None

    Raises:
    -------
    FileNotFoundError
        If the input file 'maginteract.ini' is not found.
    RuntimeError
        If previous iterations are required but not found in the HDF5 file.

    Notes:
    ------
    - Reads parameters from 'maginteract.ini' and broadcasts them across MPI processes.
    - Handles Kanamori or direct interaction Hamiltonians.
    - Checks for and uses previous iteration data if available.
    - Performs the DOS calculation for each shell and computes projected DOS for orbital contributions.

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
    basis_par = mpi.bcast(basis_par)
    magint_par = mpi.bcast(magint_par)

    mpi.report('-' * 40)
    mpi.report('Starting Density of States calculation ... ')
    mpi.report('-' * 40)

    broadening = 0.05

    # Init the SumK class
    SK = SumkDFTTools(hdf_file=general_par['filename'] + '.h5', use_dft_blocks=False)
    gf_struct = SK.gf_struct_solver_list

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
            # TODO Extend to VASP convention
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
    mpi.report('-' * 40)
    for i_sh in range(SK.n_inequiv_shells):

        # Initialize Solver
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

    mpi.report('\nChemical potential = %11.6f' % SK.chemical_potential)
    mpi.report('-' * 40)

    SK.set_Sigma([S[i_sh].Sigma_iw for i_sh in range(SK.n_inequiv_shells)])

    # Density:
    for i_sh in range(SK.n_inequiv_shells):
        S[i_sh].G_iw <<= SK.extract_G_loc()[i_sh]

        # Calculate effective atomic levels:
        eal = SK.eff_atomic_levels()
        for s in eal[i_sh]:
            eal[i_sh][s].real[abs(eal[i_sh][s].real) < solver_par['tol_eal']] = 0
            eal[i_sh][s].imag[abs(eal[i_sh][s].imag) < solver_par['tol_eal']] = 0

        # eal[i_sh] = read_eal(general_par['folder'] + '/eal_last_site_' + str(i_sh) + '.dat', nlms[i_sh])
        S[i_sh].set_atomic_levels(eal=eal[i_sh])

        mpi.report('Calculating the DOS for shell ' + str(i_sh))
        # solve it:
        if solver_par['kanamori'] is False:
            S[i_sh].GF_realomega(ommin=solver_par['ommin'], ommax=solver_par['ommax'], N_om=2 * solver_par['n_omega'],
                                 U_int=solver_par['u_int'][i_sh], J_hund=solver_par['j_hund'][i_sh])
        else:
            S[i_sh].GF_realomega(ommin=solver_par['ommin'], ommax=solver_par['ommax'], N_om=2 * solver_par['n_omega'],
                                 u4ind=U_full[i_sh])

        S[i_sh].Sigma_w << mpi.bcast(S[i_sh].Sigma_w)

    SK.put_Sigma(Sigma_imp=[S[i_sh].Sigma_w for i_sh in range(SK.n_inequiv_shells)])

    if general_par['dft_exec'] == 'Vasp':
        DOS, DOSproj, DOSproj_orb = SK.dos_wannier_basis(broadening=solver_par['broadening'], with_dc=True,
                                                         with_Sigma=True)
    else:
        DOS, DOSproj, DOSproj_orb = SK.dos_wannier_basis(broadening=solver_par['broadening'])
