# ------------------------------------------------------------------------------------#
# MagInt library
# Written by  : Leonid V. Pourovskii (CPHT Ecole Polytechnique) 2012-2025
#             : Dario Fiore Mosca (CPHT Ecole Polytechnique) 2023-2025
# Email: leonid@cpht.polytechnique.fr
# ------------------------------------------------------------------------------------#
#
#    MagInt implements the force-theorem Hubbard-I (FT-HI) approach to intersite
#    exchange interactions in correlated insulators.
#    The formalism is given in L. V. Pourovskii Phys. Rev. B 94, 115117 (2016)
#    This Python-3 version is based on TRIQS library
#
# ------------------------------------------------------------------------------------#
from triqs.operators.util.U_matrix import *
from triqs.gf import *
from hubbardI.hubbard_I import gf_hi_fullu
import triqs.utility.mpi as mpi
import numpy as np
import functools
import copy


class Solver:
    """
       Hubbard I Solver
    """

    # initialisation:
    def __init__(self, beta, l, n_iomega=1025, use_spin_orbit=False, Nmoments=5):
        """
        Initialize the solver.

        Parameters
        ----------
        beta : float
            Inverse temperature.
        l : int
            Angular momentum quantum number
        n_iomega : int, optional
            Number of Matsubara frequencies used for the Green's functions. Default is 1025.
        use_spin_orbit : bool, optional
            Whether Spin-Orbit coupling is included. Default is False.
        Nmoments : int, optional
            Number of moments. Default is 5.

        """

        self.name = "Fortran - Hubbard I Solver"
        self.beta = beta
        self.n_iomega = n_iomega
        self.UseSpinOrbit = use_spin_orbit
        self.Converged = False
        self.n_spin = 2
        self.Nmoments = Nmoments
        # self.gf_struct = fix_gf_struct_type(gf_struct)
        self.l = l
        self.Nlm = 2 * l + 1
        if (use_spin_orbit):
            # no blocks!
            self.gf_struct = [('ud', list(range(2 * self.Nlm)))]
        else:
            # up/down blocks:
            self.gf_struct = [('up', list(range(self.Nlm))), ('down', list(range(self.Nlm)))]

        g_iw_list = []

        name_list = [block for block, block_size in self.gf_struct]
        for block, block_size in self.gf_struct:
            if (use_spin_orbit):
                g_iw_list.append(
                    GfImFreq(beta=self.beta, n_points=self.n_iomega, target_shape=(len(block_size), len(block_size))))
            else:
                g_iw_list.append(
                    GfImFreq(beta=self.beta, n_points=self.n_iomega, target_shape=(len(block_size), len(block_size))))

        self.G_iw = BlockGf(name_list=name_list, block_list=g_iw_list)

        # construct Greens functions:
        self.a_list = [a for a, al in self.gf_struct]
        self.G_iw_Old = self.G_iw.copy()
        self.G0_iw = self.G_iw.copy()
        self.Sigma_iw = self.G_iw.copy()
        self.Sigma_iw_Old = self.G_iw.copy()
        self.tail = make_zero_tail(self.G_iw, self.Nmoments)

        # prepare self.ealmat
        self.ealmat = np.zeros([self.Nlm * self.n_spin, self.Nlm * self.n_spin], np.complex_)

        # Define Atomic Levels Dictionary according to the GF Bloc Structure
        self.Eff_Atomic_Levels = {}
        for a, al in self.gf_struct:
            if (self.UseSpinOrbit):
                self.Eff_Atomic_Levels[a] = np.zeros([self.Nlm * 2, self.Nlm * 2], np.complex_)
            else:
                self.Eff_Atomic_Levels[a] = np.zeros([self.Nlm, self.Nlm], np.complex_)

    def solve(self, U_int=None, J_hund=None, T=None, verbosity=0, Iteration_Number=1, Test_Convergence=0.0001, n_lev=0,
              remove_split=False, u4ind=None):
        """
        Calculate the impurity Green's function using the Hubbard-I approximation.

        Parameters
        ----------
        U_int : float, optional
            Interaction strength U.
        J_hund : float, optional
            Hund's coupling J.
        T : float, optional
            Temperature of the system.
        verbosity : int, optional
            Level of verbosity for the solver. Default is 0.
        Iteration_Number : int, optional
            Current iteration number. Default is 1.
        Test_Convergence : float, optional
            Convergence criterion for the self-consistency loop. Default is 0.0001.
        n_lev : int, optional
            Default is 0.
        remove_split : bool, optional
            Default is False.
        u4ind : ndarray, optional
            Full 4-index interaction tensor U. If provided, U_int and J_hund should be None.

        Returns
        -------
        None
            Modifies internal state of the object but does not return a value.

        Raises
        ------
        AssertionError
            If neither (U_int, J_hund) nor the full U-matrix u4ind are provided.

        Notes
        -----
        The function sets the Green's function, self-energy, etc., of the class instance.
        """

        if self.Converged:
            mpi.report("Solver %(name)s has already converged: SKIPPING" % self.__dict__)
            return

        if not ((isinstance(U_int, float) and isinstance(J_hund, float)) or isinstance(u4ind, np.ndarray)):
            mpi.report(U_int, J_hund)
            mpi.report("You need to provide either (U_int, J_hund) or the full U-matrix u4ind to %s!!" % (self.name))
            assert 0

        if mpi.is_master_node():
            self.verbosity = verbosity
        else:
            self.verbosity = 0

        nlm = self.Nlm
        if (self.UseSpinOrbit):
            nlmtot = nlm * 2
            nlms = nlmtot
        else:
            nlmtot = nlm
            nlms = nlmtot * 2

        if isinstance(u4ind, np.ndarray):
            # set 4-ind and 2-ind U from input
            U = u4ind
            ummss = np.zeros((nlms, nlms))
            for i in range(nlms):
                for j in range(nlms):
                    ummss[i, j] = U[i, j, i, j].real - U[i, j, j, i].real
        else:
            # construct U from input U_int and J_hund
            ur, umn, ujmn = self.__set_umatrix(U=U_int, J=J_hund, T=T)
            ummss = np.zeros((nlms, nlms))
            U = np.zeros((nlms, nlms, nlms, nlms))
            ummss[0:nlm, 0:nlm] = ujmn
            ummss[nlm:nlms, nlm:nlms] = ujmn
            ummss[0:nlm, nlm:nlms] = umn
            ummss[nlm:nlms, 0:nlm] = umn
            for is1 in range(2):
                for is2 in range(2):
                    U[is1 * nlm:(is1 + 1) * nlm, is2 * nlm:(is2 + 1) * nlm, is1 * nlm:(is1 + 1) * nlm,
                    is2 * nlm:(is2 + 1) * nlm] = ur

        self.ummss = ummss

        M = [x for x in self.G_iw.mesh]
        self.zmsb = np.array([x for x in M], np.complex_)

        self.__save_eal('eal.dat', Iteration_Number)

        mpi.report("Starting Fortran solver %(name)s" % self.__dict__)

        self.Sigma_iw_Old <<= self.Sigma_iw
        self.G_iw_Old <<= self.G_iw

        iwmax = int(len(self.G_iw.mesh))
        gf = np.zeros((nlm * self.n_spin, nlm * self.n_spin, iwmax), dtype=np.cdouble, order='F')
        tail = np.zeros((self.Nmoments, nlm * self.n_spin, nlm * self.n_spin), dtype=np.cdouble, order='F')

        # call the fortran solver:
        temp = 1.0 / self.beta

        self.atocc, self.atmag = gf_hi_fullu(gf=gf, tail=tail, e0f=self.ealmat, u=U, ummss=ummss, zmsb=self.zmsb,
                                             nlm=nlm, iwmax=iwmax, nmom=self.Nmoments, ns=self.n_spin, temp=temp,
                                             verbosity=self.verbosity, remove_split=remove_split, n_lev=n_lev)

        # self.sig = sigma_atomic_fullu(gf=self.gf,e0f=self.eal,zmsb=self.zmsb,ns=self.n_spin,nlm=self.Nlm)
        if (self.verbosity == 0):
            # No fortran output, so give basic results here
            mpi.report("Atomic occupancy in Hubbard I Solver  : %s" % self.atocc)
            mpi.report("Atomic magn. mom. in Hubbard I Solver : %s" % self.atmag)

        # transfer the data to the GF class:
        M = {}
        isp = -1
        for a, al in self.gf_struct:
            isp += 1
            M[a] = np.array(gf[isp * nlmtot:(isp + 1) * nlmtot, isp * nlmtot:(isp + 1) * nlmtot, :]).transpose(2, 0,
                                                                                                               1).copy()
            # Tails
            for i in range(self.Nmoments - 1):
                self.tail[isp][i + 1] = tail[i][isp * nlmtot:(isp + 1) * nlmtot, isp * nlmtot:(isp + 1) * nlmtot]

        # glist = lambda : [ GfImFreq(indices = al, beta = self.beta, n_points = self.n_iomega) for a,al in self.gf_struct]
        # self.G = BlockGf(name_list = self.a_list, block_list = glist(),make_copies=False)
        glist = []
        name_list = [block for block, block_size in self.gf_struct]
        for block, block_size in self.gf_struct:
            glist.append(
                GfImFreq(beta=self.beta, n_points=self.n_iomega, target_shape=(len(block_size), len(block_size))))

        self.G_iw = BlockGf(name_list=name_list, block_list=glist)
        self.__copy_Gf(self.G_iw, M)

        # Self energy:
        self.G0_iw <<= iOmega_n

        M = [self.ealmat[isp * nlmtot:(isp + 1) * nlmtot, isp * nlmtot:(isp + 1) * nlmtot] for isp in
             range((2 * self.Nlm) // nlmtot)]
        self.G0_iw -= M
        self.Sigma_iw <<= self.G0_iw - inverse(self.G_iw)

        # invert G0
        self.G0_iw.invert()

        def test_distance(G1, G2, dist):
            def f(G1, G2):
                # print abs(G1.data - G2.data)
                dS = max(abs(G1.data - G2.data).flatten())
                aS = max(abs(G1.data).flatten())
                return dS <= aS * dist

            return functools.reduce(lambda x, y: x and y, [f(g1, g2) for (i1, g1), (i2, g2) in zip(G1, G2)])

        mpi.report("\nChecking Sigma for convergence...\nUsing tolerance %s" % Test_Convergence)
        self.Converged = test_distance(self.Sigma_iw, self.Sigma_iw_Old, Test_Convergence)

        if self.Converged:
            mpi.report("Solver HAS CONVERGED")
        else:
            mpi.report("Solver has not yet converged")

    def GF_realomega(self, ommin, ommax, N_om, U_int=None, J_hund=None, T=None, verbosity=0, broadening=0.01, n_lev=0,
                     remove_split=False, u4ind=None):
        """
        Calculates the Green's Function (GF) and spectral function on the real frequency axis.

        Parameters:
        -----------
        ommin : float
            Minimum frequency for the calculation.
        ommax : float
            Maximum frequency for the calculation.
        N_om : int
            Number of frequency points.
        U_int : float, optional
            Interaction parameter U for defining the interactions in the system.
            Either (U_int, J_hund) should be provided or the full U-matrix u4ind.
        J_hund : float, optional
            Hund's coupling parameter J.
            Either (U_int, J_hund) should be provided or the full U-matrix u4ind.
        T : float, optional
            Temperature of the system. Default is None.
        verbosity : int, optional
            Level of verbosity for logging. Default is 0.
        broadening : float, optional
            Broadening factor for the GF. Default is 0.01.
        n_lev : int, optional
            Number of levels included. Default is 0.
        remove_split : bool, optional
            If True, removes the splitting in the GF. Default is False.
        u4ind : numpy.ndarray, optional
            A full 4-index U-matrix. If provided, (U_int, J_hund) are not used.

        Returns:
        --------
        None. Modifies internal class attributes to store computed GF, spectral functions, etc.

        Note:
        -----
        Either the interaction parameters (U_int, J_hund) should be provided, or the full U-matrix u4ind.
        """

        if not ((isinstance(U_int, float) and isinstance(J_hund, float)) or isinstance(u4ind, np.ndarray)):
            print(U_int, J_hund)
            mpi.report("You need to provide either (U_int, J_hund) or the full U-matrix u4ind to %s!!" % (self.name))
            assert 0

        delta_om = (ommax - ommin) / (1.0 * (N_om - 1))
        omega = np.zeros([N_om], np.complex_)

        nlm = self.Nlm
        #if (self.UseSpinOrbit):
        nlmtot = nlm * 2
        nlms = nlmtot
        #else:
        #    nlmtot = nlm
        #    nlms = nlmtot

        if isinstance(u4ind, np.ndarray):
            # set 4-ind and 2-ind U from input
            U = u4ind
            ummss = np.zeros((nlms, nlms))
            for i in range(nlms):
                for j in range(nlms):
                    ummss[i, j] = U[i, j, i, j].real - U[i, j, j, i].real
        else:
            # construct U from input U_int and J_hund
            ur, umn, ujmn = self.__set_umatrix(U=U_int, J=J_hund, T=T)
            ummss = np.zeros((nlms, nlms))
            U = np.zeros((nlms, nlms, nlms, nlms))
            ummss[0:nlm, 0:nlm] = ujmn
            ummss[nlm:nlms, nlm:nlms] = ujmn
            ummss[0:nlm, nlm:nlms] = umn
            ummss[nlm:nlms, 0:nlm] = umn
            for is1 in range(2):
                for is2 in range(2):
                    U[is1 * nlm:(is1 + 1) * nlm, is2 * nlm:(is2 + 1) * nlm, is1 * nlm:(is1 + 1) * nlm,
                    is2 * nlm:(is2 + 1) * nlm] = ur

        self.ummss = ummss

        for i in range(N_om):
            omega[i] = ommin + delta_om * i + 1j * broadening

        iwmax = int(len(self.G_iw.mesh))
        gf = np.zeros((nlm * self.n_spin, nlm * self.n_spin, iwmax), dtype=np.cdouble, order='F')
        tail = np.zeros((self.Nmoments, nlm * self.n_spin, nlm * self.n_spin), dtype=np.cdouble, order='F')

        temp = 1.0 / self.beta
        self.atocc, self.atmag = gf_hi_fullu(gf=gf, tail=tail, e0f=self.ealmat, u=U, ummss=ummss,
                                             zmsb=omega, nlm=nlm, iwmax=iwmax, nmom=self.Nmoments, ns=self.n_spin,
                                             temp=temp, verbosity=verbosity, remove_split=remove_split,
                                             n_lev=n_lev)

        # transfer the data to the GF class:
        if (self.UseSpinOrbit):
            nlmtot = self.Nlm * 2  # only one block in this case!
        else:
            nlmtot = self.Nlm

        M = {}
        isp = -1
        for a, al in self.gf_struct:
            isp += 1
            M[a] = np.array(gf[isp * nlmtot:(isp + 1) * nlmtot, isp * nlmtot:(isp + 1) * nlmtot, :]).transpose(2, 0,
                                                                                                               1).copy()
            for i in range(self.Nmoments - 1):
                self.tail[isp][i + 1] = tail[i][isp * nlmtot:(isp + 1) * nlmtot, isp * nlmtot:(isp + 1) * nlmtot]

        g_w_list = []

        name_list = [block for block, block_size in self.gf_struct]
        for block, block_size in self.gf_struct:
            if (self.UseSpinOrbit):
                self.Nlm = int(len(block_size) / 2)
                g_w_list.append(
                    GfReFreq(window=(ommin, ommax), n_points=N_om, target_shape=(len(block_size), len(block_size))))
            else:
                self.Nlm = len(block_size)
                g_w_list.append(
                    GfReFreq(window=(ommin, ommax), n_points=N_om, target_shape=(len(block_size), len(block_size))))

        self.G_w = BlockGf(name_list=name_list, block_list=g_w_list)

        self.__copy_Gf(self.G_w, M)

        # Self energy:
        self.G0_w = self.G_w.copy()
        self.Sigma_w = self.G_w.copy()
        self.G0_w <<= Omega + 1j * broadening

        M = [self.ealmat[isp * nlmtot:(isp + 1) * nlmtot, isp * nlmtot:(isp + 1) * nlmtot] for isp in
             range((2 * self.Nlm) // nlmtot)]
        self.G0_w -= M
        self.Sigma_w <<= self.G0_w - inverse(self.G_w)
        self.Sigma_w.note = 'ReFreq'  # This is important for the put_Sigma routine!!!

        # sigmamat = sigma_atomic_fullu(gf=gf,e0f=self.ealmat,zmsb=omega,nlm=self.Nlm,ns=self.n_spin)

        # return omega,gf,sigmamat

    def __save_eal(self, Filename, it):
        if mpi.is_master_node():
            f = open(Filename, 'a')
            f.write('\neff. atomic levels, Iteration %s\n' % it)
            for i in range(self.Nlm * self.n_spin):
                for j in range(self.Nlm * self.n_spin):
                    f.write("%10.6f %10.6f   " % (self.ealmat[i, j].real, self.ealmat[i, j].imag))
                f.write("\n")
            f.close()

    def __copy_Gf(self, G, data):
        """ Copies data and tail to Gf object GF """
        for s, g in G:
            g.data[:, :, :] = data[s][:, :, :]
        #   for imom in range(1,min(self.Nmoments,8)):
        #       g.tail.data[1+imom,:,:]=tail[s][imom]

    def set_atomic_levels(self, eal):
        """
        Sets atomic level data based on the provided dictionary.

        This method updates the atomic level matrix (self.ealmat) and
        effective atomic levels (self.Eff_Atomic_Levels) attributes
        using the dictionary of atomic levels provided.

        Parameters:
        -----------
        eal : dict
            Dictionary containing atomic levels data.
            Keys are indices, and values are 2D arrays or matrices
            representing the atomic level information for the respective index.

        Returns:
        --------
        None. Modifies internal class attributes `self.ealmat` and `self.Eff_Atomic_Levels`.
        """

        # Function code ...

        assert (type(eal) == dict), "Give a dictionary to set_atomic_levels!"

        cnt = 0
        self.ealmat[:, :] *= 0.0

        for ind in eal:
            self.Eff_Atomic_Levels[ind] = copy.deepcopy(eal[ind])

            if self.UseSpinOrbit:
                for ii in range(self.Nlm * 2):
                    for jj in range(self.Nlm * 2):
                        self.ealmat[ii, jj] = self.Eff_Atomic_Levels[ind][ii, jj]
            else:
                for ii in range(self.Nlm):
                    for jj in range(self.Nlm):
                        self.ealmat[cnt * self.Nlm + ii, cnt * self.Nlm + jj] = self.Eff_Atomic_Levels[ind][ii, jj]

            cnt += 1

    def __set_umatrix(self, U, J, T=None):
        # U matrix:
        # l = (Nlm-1)/2
        # If T is specified, it is used to transform the Basis set
        if self.l > 0:
            Umat = U_matrix(l=self.l, U_int=U, J_hund=J, basis='spherical', T=T)
            U, Up = reduce_4index_to_2index(Umat)
        else:
            Umat = np.zeros((1, 1, 1, 1), dtype=float)
            Umat[0, 0, 0, 0] = U
            U = np.zeros((1, 1), dtype=float)
            Up = np.zeros((1, 1), dtype=float)
            Up[0, 0] = U

        return Umat, Up, U
