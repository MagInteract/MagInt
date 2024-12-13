from triqs_dft_tools.sumk_dft import *
from triqs_dft_tools.converters.wien2k import *
from MagInt.HubbardI_solver import Solver as Solver
from utils import *
import triqs.utility.mpi as mpi
from h5 import *
import numpy as np
import os


dft_filename = os.getcwd().rpartition('/')[2]
dft_exec = 'Wien2k'

# Parameters for the Hubbard U initialization
U_int = [3.2, 3.2]  # Interaction strength
J_hund = [0.5, 0.5]  # Hund's coupling
Kanamori = False  # If True, use Kanamori parameterization for U

# Solver Parameters
beta = 40
Loops = 7  # Number of DMFT sc-loops
mixing = 0.7  # Mixing factor
DC_type = 0  # 0...FLL, 1...Held, 2... AMF, 3...Lichtenstein
chemical_potential_init = 0.0  # initial chemical potential
split = 0.0000005
verbosity = 2
use_spin_orbit = True
n_mult = [1, 1]
n_ineq = 0
LS = [2, 2]  # Orbital momentum for solver for each shell
n_el = [1, 1]  # Number of electrons

# Convert DMFT input:
Converter = Wien2kConverter(filename=dft_filename)
Converter.convert_dft_input()
mpi.barrier()

# Init the SumK class
SK = SumkDFT(hdf_file=dft_filename + '.h5', use_dft_blocks=False)
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
    nlm.append(2 * LS[i_sh] + 1)
    nlms.append(2 * nlm[i_sh])
    if Kanamori is True:
        # TODO Extend to VASP convention
        U_full_val, U_dc_val, J_dc_val = Kanamori_hamiltonian(U_int=U_int[i_sh], J_hund=J_hund[i_sh], orbital=LS[i_sh],
                                                              dft_exec=dft_exec)
        U_full.append(U_full_val)
        mpi.report('Generating the Kanamori U from U_int : ', U_int[i_sh], ' and J_hund : ', J_hund[i_sh])
    else:
        U_dc_val = U_int[i_sh]
        J_dc_val = J_hund[i_sh]
        mpi.report('Generating the full U from U_int : ', U_int[i_sh], ' and J_hund : ', J_hund[i_sh])

    U_dc.append(U_dc_val)
    J_dc.append(J_dc_val)

    # Initialization DC
    useval.append(U_dc[i_sh] * (n_el[i_sh] - 0.5) - J_dc[i_sh] * (n_el[i_sh] * 0.5 - 0.5))

# check if there are previous runs:
previous_runs = 0
previous_present = False
if mpi.is_master_node():
    f = HDFArchive(dft_filename + '.h5', 'a')
    if 'dmft_output' in f:
        ar = f['dmft_output']
        if 'iterations' in ar:
            previous_present = True
            previous_runs = ar['iterations']
    else:
        f.create_group('dmft_output')
    del f
previous_runs = mpi.bcast(previous_runs)
previous_present = mpi.bcast(previous_present)

# Setup Solver
mpi.report('Sumk to Solver: %s' % SK.sumk_to_solver)
mpi.report('GF struct sumk: %s' % SK.gf_struct_sumk)
mpi.report('GF struct solver: %s' % SK.gf_struct_solver)

chemical_potential = chemical_potential_init

# Init the Hubbard-I solver:
S = []
dm = []
for i_sh in range(SK.n_inequiv_shells):
    S.append(Solver(beta=beta, l=LS[i_sh], use_spin_orbit=True))

    iteration_offset = 0

    # load previous data: old self-energy, chemical potential, DC correction
    if previous_present:
        mpi.report('PREVIOUS PRESENT')
        if mpi.is_master_node():
            ar = HDFArchive(dft_filename + '.h5', 'a')
            S[i_sh].Sigma_iw << ar['dmft_output']['Sigma_iw_%s' % i_sh]
            del ar
            SK.chemical_potential, SK.dc_imp, SK.dc_energ = SK.load(['chemical_potential', 'dc_imp', 'dc_energ'])
        S[i_sh].Sigma_iw << mpi.bcast(S[i_sh].Sigma_iw)
        SK.chemical_potential = mpi.bcast(SK.chemical_potential)
        SK.dc_imp = mpi.bcast(SK.dc_imp)
        SK.dc_energ = mpi.bcast(SK.dc_energ)

    # calculated fixed DC at the beginning
    dm.append(S[i_sh].G_iw.density())

    SK.calc_dc(dm[i_sh], U_interact=U_int[i_sh], J_hund=J_hund[i_sh], orb=i_sh, use_dc_formula=DC_type,
               use_dc_value=useval[i_sh])
    if not previous_present:
        S[i_sh].Sigma_iw << SK.dc_imp[i_sh]['ud'][0, 0]

    mpi.report('%s DMFT cycles requested. Starting with iteration %s.' % (Loops, iteration_offset))

# DMFT loop:
for iteration_number in range(1, Loops + 1):

    itn = iteration_number + previous_runs

    # put Sigma_iw into the SumK class:
    SK.set_Sigma([S[i].Sigma_iw for i in range(SK.n_inequiv_shells)])

    # Compute the SumK, possibly fixing mu by dichotomy
    chemical_potential = SK.calc_mu(precision=0.000001)

    for i_sh in SK.corr_to_inequiv:

        # Density:
        S[i_sh].G_iw <<= SK.extract_G_loc()[i_sh]
        mpi.report("Total charge of Gloc : ", S[i_sh].G_iw.total_density().real)

        dm = S[i_sh].G_iw.density()
        mpi.report("Orbital densities of local Green function of site %s:" % (i_sh))
        for s in dm:
            mpi.report("Block %s: " % s)
            for ii in range(len(dm[s])):
                st = ''
                for jj in range(len(dm[s])):
                    if (dm[s][ii, jj].real > 0):
                        st += "   %.4f" % (dm[s][ii, jj].real)
                    else:
                        st += "  %.4f" % (dm[s][ii, jj].real)
                mpi.report(st)

        eal = SK.eff_atomic_levels()[i_sh]

        tol_eal = 0.0005
        for s in eal:
            eal[s].real[abs(eal[s].real) < tol_eal] = 0
            eal[s].imag[abs(eal[s].imag) < tol_eal] = 0

        for s in eal:
            mpi.report("\nEffective atomic levels of site %s:" % (i_sh))
            mpi.report("Block real %s:" % (i_sh))
            for ii in range(len(eal[s])):
                st = ''
                for jj in range(len(eal[s])):
                    if eal[s][ii, jj].real > 0:
                        st += "   %.4f" % (eal[s][ii, jj].real)
                    else:
                        st += "  %.4f" % (eal[s][ii, jj].real)
                mpi.report(st)
            mpi.report("Block Imag %s:" % (i_sh))
            for ii in range(len(eal[s])):
                st = ''
                for jj in range(len(eal[s])):
                    if eal[s][ii, jj].real > 0:
                        st += "   %.4f" % (eal[s][ii, jj].imag)
                    else:
                        st += "  %.4f" % (eal[s][ii, jj].imag)
                mpi.report(st)

        S[i_sh].set_atomic_levels(eal=eal)

        # solve it:
        if mpi.is_master_node():
            if os.path.exists('eal_site_%s.dat' % i_sh) and itn > 1:
                os.rename('eal_site_%s.dat' % i_sh, 'eal.dat')
            else:
                if os.path.exists('eal.dat'): os.remove('eal.dat')

        # solve it:
        if Kanamori is True:
            S[i_sh].solve(u4ind=U_full[i_sh], Test_Convergence=1e-30, verbosity=verbosity)
        else:
            S[i_sh].solve(U_int=U_int[i_sh], J_hund=J_hund[i_sh], Test_Convergence=1e-30, verbosity=verbosity)

        mpi.barrier()

        if mpi.is_master_node():
            if os.path.exists('eal.dat'):
                os.rename('eal.dat', 'eal_site_%s.dat' % i_sh)
            if os.path.exists('atomic_levels.dat'):
                os.rename('atomic_levels.dat', 'atomic_levels_site_%s.dat' % i_sh)
            if os.path.exists('states.dat'):
                os.rename('states.dat', 'states_site_%s.dat' % i_sh)

    # Now mix Sigma_iw and G with factor Mix, if wanted:
    if itn > 1 or previous_present:
        if mpi.is_master_node() and (mixing < 1.0):
            ar = HDFArchive(dft_filename + '.h5', 'a')
            mpi.report("Mixing Sigma_iw and G with factor %s" % mixing)
            for i_sh in SK.corr_to_inequiv:
                S[i_sh].Sigma_iw << mixing * S[i_sh].Sigma_iw + (1.0 - mixing) * ar['dmft_output']['Sigma_iw_%s' % i_sh]
                S[i_sh].G_iw << mixing * S[i_sh].G_iw + (1.0 - mixing) * ar['dmft_output']['G_iw_%s' % i_sh]
            del ar
        for i_sh in SK.corr_to_inequiv:
            S[i_sh].G_iw << mpi.bcast(S[i_sh].G_iw)
            S[i_sh].Sigma_iw << mpi.bcast(S[i_sh].Sigma_iw)

    # correlation energy calculations:
    SK.correnerg = 0.0
    for i_sh in SK.corr_to_inequiv:
        SK.correnerg += 0.5 * (S[i_sh].G_iw * S[i_sh].Sigma_iw).total_density() * n_mult[i_sh]
    mpi.report("Correlation energy = %s" % SK.correnerg.real)

    # store the impurity self-energy, GF as well as correlation energy in h5
    if mpi.is_master_node():
        ar = HDFArchive(dft_filename + '.h5', 'a')
        ar['dmft_output']['iterations'] = itn
        for i_sh in SK.corr_to_inequiv:
            ar['dmft_output']['G_iw_%s' % i_sh] = S[i_sh].G_iw
            ar['dmft_output']['Sigma_iw_%s' % i_sh] = S[i_sh].Sigma_iw
        del ar

    # Save essential SumkDFT data:
    SK.save(['chemical_potential', 'dc_imp', 'dc_energ', 'correnerg'])
    if mpi.is_master_node():
        for i_sh in SK.corr_to_inequiv:
            print('DC after solver for shell', i_sh, ' : ', SK.dc_imp[i_sh])
            mpi.report("Total charge of impurity problem : %.6f" % S[i_sh].G_iw.total_density())

# find exact chemical potential
SK.chemical_potential = SK.calc_mu(precision=0.000001)

# calculate and save occupancy matrix in the Bloch basis for Wien2k charge denity recalculation
dN, d = SK.calc_density_correction(filename=dft_filename + '.qdmft')

mpi.report("Trace of Density Matrix: %s" % d)

# store correlation energy contribution to be read by Wien2ki and then included to DFT+DMFT total energy
if (mpi.is_master_node()):
    SK.correnerg -= SK.dc_energ[0]
    f = open(dft_filename + '.qdmftup', 'a')
    f.write("%.16f\n" % (SK.correnerg * 0.5))
    f.close()
    f = open(dft_filename + '.qdmftdn', 'a')
    f.write("%.16f\n" % (SK.correnerg * 0.5))
    f.close()

