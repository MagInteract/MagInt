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
from triqs_dft_tools.sumk_dft import SumkDFT
import triqs.utility.mpi as mpi
from MagInt.HubbardI_interact import *
from MagInt.Multipolar import *
from MagInt.MagInteract import *
from MagInt.Read_input import *
from MagInt.utils import *
from MagInt import *
import numpy as np
import sys


def run_magint():
    """
    Performs the MagInt Intersite Exchange Interactions (IEI) calculation.

    This function initializes parameters, sets up solvers, handles previous iterations if available,
    and calculates magnetic interactions using input from a configuration file and interaction parameters.

    Parameters:
    -----------
    None

    Raises:
    -------
    FileNotFoundError
        If the input file 'maginteract.ini' is not found.
    RuntimeError
        If previous iterations are required but not found in the HDF5 file.
    AssertionError
        If the basis parameters are not correctly formatted or have unsupported types.

    Notes:
    ------
    - Reads parameters from 'maginteract.ini' and broadcasts them across MPI processes.
    - Sets up the `SumkDFT` object and initializes solvers for each inequivalent shell.
    - Handles Kanamori or direct interaction Hamiltonians.
    - Checks for and uses previous iteration data if available.
    - Calculates effective atomic levels and initializes Green's functions.
    - Computes and transforms magnetic interaction matrices.
    - Handles both simple J and combined J operator basis for IEI calculations.

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
    mpi.report('Starting MagInteract calculation ... ')
    mpi.report('-' * 40)

    # Init the SumK class
    SK = SumkDFT(hdf_file=general_par['filename'] + '.h5', use_dft_blocks=False)
    if not general_par['use_symmetries']: SK.symm_op = 0

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

    # Magnetic interactions to be calculated
    label_corrsite = {}
    GSM = {}
    StBas_item = {}
    S_HI = {}

    for i_sh in range(SK.n_inequiv_shells):

        if isinstance(magint_par['included_atoms'], list) and i_sh not in magint_par['included_atoms']: continue

        label_corrsite[magint_par['correlated_atoms'][i_sh]] = i_sh
        if type(basis_par['ps_j'][i_sh]) is list:
            assert len(basis_par['ps_j'][i_sh]) == 2, \
                'Value of ps_j for site %s is %s. ps_j for a site can be either a float or a list of two floats' % (
                    ish, basis_par['ps_j'][i_sh])
            if basis_par['action'] == 'merge':
                GSM[magint_par['correlated_atoms'][i_sh]] = int(
                    2 * basis_par['ps_j'][i_sh][0] + 2 * basis_par['ps_j'][i_sh][1] + 2)
            elif basis_par['action'] == 'dirprod':
                GSM[magint_par['correlated_atoms'][i_sh]] = int(
                    (2 * basis_par['ps_j'][i_sh][0] + 1) * (2 * basis_par['ps_j'][i_sh][1] + 1))
        else:
            GSM[magint_par['correlated_atoms'][i_sh]] = int((2 * basis_par['ps_j'][i_sh]) + 1)
        StBas_item[magint_par['correlated_atoms'][i_sh]] = ['Ion' + str(i_sh)]

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

    # Calculate effective Atomic levels
    eal = SK.eff_atomic_levels()

    ifSO=solver_par['use_spin_orbit']

    # Set DC
    mpi.report('-' * 40)
    mpi.report('DC after reading SK: ')
    for i_sh in range(SK.n_inequiv_shells):
        if ifSO: 
            SK.dc_imp[i_sh]['ud'][0, 0] = useval[i_sh]
        else:
            SK.dc_imp[i_sh]['up'][0, 0] = useval[i_sh]
            SK.dc_imp[i_sh]['down'][0, 0] = useval[i_sh]
        if mpi.is_master_node() and general_par['verbosity'] > 1:
            if ifSO:
                print('DC after reading SK for shell :', i_sh, '\n', SK.dc_imp[i_sh]['ud'])
            else:
                print('DC after reading SK for shell :', i_sh, '\n', SK.dc_imp[i_sh]['up'])

    # Init the Solver:
    S = []
    dm = {}
    if ifSO:
        dm['ud'] = np.zeros((nlms[0], nlms[0]), complex)
    else:
        dm['up'] = np.zeros((nlm[0], nlm[0]), complex)
        dm['down'] = np.zeros((nlm[0], nlm[0]), complex)
    for i_sh in range(SK.n_inequiv_shells):
        # Initialize Solver
        S.append(Solver(beta=solver_par['beta'], l=solver_par['ls'][i_sh], n_iomega=solver_par['n_iomega'],
                        use_spin_orbit=solver_par['use_spin_orbit']))
        SK.calc_dc(dm, U_interact=solver_par['u_int'][i_sh], J_hund=solver_par['j_hund'][i_sh], orb=i_sh,
                   use_dc_formula=solver_par['dc_type'], use_dc_value=useval[i_sh])

        if not previous_present:
            if ifSO:
                S[i_sh].Sigma_iw << SK.dc_imp[0]['ud'][0, 0]
            else:
                S[i_sh].Sigma_iw << SK.dc_imp[0]['up'][0, 0]

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

        # Calculate effective atomic levels:
        eal = SK.eff_atomic_levels()
        for s in eal[i_sh]:
            eal[i_sh][s].real[abs(eal[i_sh][s].real) < solver_par['tol_eal']] = 0
            eal[i_sh][s].imag[abs(eal[i_sh][s].imag) < solver_par['tol_eal']] = 0

        #if mpi.is_master_node():
        #    print_arr(eal[i_sh]['ud'], log='Effective Atomic Levels : ', decd=6,
        #              fname='eal_in_mag_int_' + str(i_sh) + '.dat',
        #              prn_zero_imag=False)
        #    print_arr(eal[i_sh]['ud'].real, log='Effective Atomic Levels : ', decd=6, prn_zero_imag=False)

        # Exclude non-interacting ions from the Solver
        if isinstance(magint_par['included_atoms'], list) and i_sh not in magint_par['included_atoms']: continue
        # Initialise the Hubbard-I solver for each correlated ion
        S_HI[magint_par['correlated_atoms'][i_sh]] = S[i_sh]

    SK.set_Sigma([S[i_sh].Sigma_iw for i_sh in range(SK.n_inequiv_shells)])
    # Initialize MagInt
    MagInt = MagInteract(general_par, magint_par, label_corrsite, GSM, SK, S_HI)

    # initialize HubbardI_interact "solver"
    S_INT = {}
    for num, at_type in enumerate(MagInt.interact_types):
        StBas = MagInt.Retrieve_StBas(StBas_item[at_type], h5name=general_par['folder'] + '/Standard_Basis.h5')
        S_INT[at_type] = []
        for val in MagInt.interact_sites[at_type]:
            icrsh = val[0]
            ieq = val[1]
            ind0 = SK.inequiv_to_corr[icrsh]
            ish_SK = ind0 + ieq
            nlm_sh = 2 * solver_par['ls'][icrsh] + 1
            mpi.report("\nSet up standard basis for type %s and equivalent site %s" % (at_type, ieq))
            StBas_site = MagInt.rot_stbas(StBas, solver_par['n_el'][icrsh], nlm_sh, SK.rot_mat[ish_SK],
                                          SK.rot_mat_time_inv[ish_SK])
            StBas_site = np.asfortranarray(StBas_site)
            if solver_par['kanamori'] is True:
                S_INT[at_type].append(
                    HubbardI_interact(beta=solver_par['beta'], l=solver_par['ls'][icrsh],
                                      n_iomega=solver_par['n_iomega'], n_lev=GSM[at_type],
                                      u4ind=U_full[icrsh], verbosity=2, st_bas=StBas_site, CalcOvl=True))
            else:
                S_INT[at_type].append(
                    HubbardI_interact(beta=solver_par['beta'], l=solver_par['ls'][icrsh],
                                      n_iomega=solver_par['n_iomega'], n_lev=GSM[at_type],
                                      U_int=solver_par['u_int'][icrsh], J_hund=solver_par['j_hund'][icrsh],
                                      verbosity=general_par['verbosity'], st_bas=StBas_site, CalcOvl=True))
                # verbosity=general_par['verbosity'], st_bas=StBas, CalcOvl=True))
            S_INT[at_type][-1].set_ud_levels(eal=eal[icrsh], rmat=SK.rot_mat[ish_SK],
                                             rmat_time_inv=SK.rot_mat_time_inv[ish_SK])

    V, V_4ind, pairs = MagInt.calc_MagInt(general_par, S_INT, calc_off_diag=True)

    mpi.report("Transformation to the standard basis")
    MagInt.trans_V(V, V_4ind, S_INT)
    mpi.report("Transformation ends")

    if mpi.is_master_node():
        for num, at_type in enumerate(MagInt.interact_types):
            for val in MagInt.interact_sites[at_type]:
                ieq = val[1]
                mpi.report('\nOverlap matrix atom %s and site %s :' % (at_type, ieq))
                ovlmat = S_INT[at_type][ieq].ovlmat
                for m in range(GSM[at_type]):
                    str_1 = ' '
                    for m1 in range(GSM[at_type]):
                        str_1 += '   %8.4f %8.4f' % (ovlmat[m, m1].real, ovlmat[m, m1].imag)
                    mpi.report(str_1)
                unitary = is_unitary(ovlmat)
                if unitary:
                    mpi.report('The Matrix is unitary')
                else:
                    mpi.report('The Matrix is not unitary')

        # index basis initialization
        try:
            at1 = magint_par['correlated_atoms'].index(magint_par['atom1'])
            at2 = magint_par['correlated_atoms'].index(magint_par['atom2'])
        except ValueError:
            print(f"Correlated atoms not found. Check the maginteract.ini file")

        if type(basis_par['ps_j'][at1]) is float and basis_par['ps_j'][at1] == basis_par['ps_j'][at2]:
            VM = Multipolar(J=basis_par['ps_j'][at1], basis='real', print_mat=True,
                            conv_mon_dipol=magint_par['mult_conv'])
        elif type(basis_par['ps_j'][at1]) is float and basis_par['ps_j'][at1] != basis_par['ps_j'][at2]:
            VM = Multipolar(J=basis_par['ps_j'][at1], basis='real', print_mat=True,
                            conv_mon_dipol=magint_par['mult_conv'], J1=basis_par['ps_j'][at2])
        elif len(basis_par['ps_j'][at1]) == 2:
            # combined operator basis
            VM1 = Multipolar(J=basis_par['ps_j'][0][0], basis='real', print_mat=False,
                             conv_mon_dipol=magint_par['mult_conv'])
            VM2 = Multipolar(J=basis_par['ps_j'][0][1], basis='real', print_mat=False,
                             conv_mon_dipol=magint_par['mult_conv'])
            #
            VM = Multipolar(Mult_Comp=[VM1, VM2], action=basis_par['action'], print_mat=True)
        else:
            assert 0, "The basis parameter 'ps_j' should be either a float or a list"

        VM.Mult_interact(V_4ind, MagInt, fname='V_Mult_%s-%s.dat' % (magint_par['atom1'], magint_par['atom2']),
                         tol_prnt=magint_par['tol_prnt'])
        MI_h5name = 'Vmult'
        VM.store_results_in_h5(MagInt, MI_h5name)
