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
from triqs.gf import *
from triqs.operators.util.U_matrix import *
from MagInt.HubbardI_solver import Solver
from MagInt.utils import *
from hubbardI.hubbard_I import gf_hi_fullu_int, ops_on_state
import triqs.utility.mpi as mpi
import numpy as np
from scipy.special import factorial

from h5 import *


class HubbardI_interact(Solver):
    """
      Version of Hubbard-I to produce Sigmas for each of atomic configuration in the GS multiplet.
      Finally calculates differences between S.Sigma Hub-I and those self-energies that are used in magnetic interaction
      calculations
    """

    def __init__(self, beta, l, n_lev, U_int=None, J_hund=None, u4ind=None, T=None, n_iomega=1025,
                 use_spin_orbit=True, lad_op=None, st_bas=None, n_bas=None, gs_occ=None, CalcOvl=False,
                 Nmoments=5, verbosity=0):
        """
        Initialize the HubbardI S instance with given parameters.

        Parameters:
        -----------
        beta : float
            Inverse temperature.
        l : int
            Angular momentum quantum number.
        n_lev : int
            Number of levels.
        U_int : float, optional
            Interaction term. Default is None.
        J_hund : float, optional
            Hund's coupling term. Default is None.
        u4ind : np.ndarray, optional
            Full U-matrix. Default is None.
        T : float, optional
            Temperature. Default is None.
        n_iomega : int, optional
            Number of frequency points. Default is 1025.
        use_spin_orbit : bool, optional
            Whether to use spin orbit coupling. Default is True.
        lad_op : type, optional
            Default is None.
        st_bas : type, optional
            Default is None.
        n_bas : type, optional
            Default is None.
        gs_occ : type, optional
            Default is None.
        CalcOvl : bool, optional
            Whether to calculate overlap. Default is False.
        Nmoments : int, optional
            Number of moments. Default is 5.
        verbosity : int, optional
            Verbosity level. Default is 0.

        Notes:
        ------
        Verbosity level will be set to 0 for non-master nodes.
        """

        super().__init__(beta, l, n_iomega=n_iomega, use_spin_orbit=use_spin_orbit, Nmoments=Nmoments)
        self.n_lev = n_lev
        self.use_spin_orbit = use_spin_orbit
        self.st_bas = st_bas
        self.gs_occ = gs_occ
        self.n_bas = n_bas
        self.lad_op = lad_op
        self.CalcOvl = CalcOvl
        # if not ((U_int and J_hund) or isinstance(u4ind,np.ndarray)):
        if not ((isinstance(U_int, float) and isinstance(J_hund, float)) or isinstance(u4ind, np.ndarray)):
            mpi.report("You need to provide either (U_int, J_hund) or the full U-matrix u4ind to HubI initialization!!")
            assert 0
        self.u4ind = u4ind
        self.l = l
        self.U_int = U_int
        self.J_hund = J_hund
        self.T = T
        if mpi.is_master_node():
            self.verbosity = verbosity
        else:
            self.verbosity = 0

        # self.gf_struct = fix_gf_struct_type(gf_struct)
        # for block, block_size in self.gf_struct:
        #     if (use_spin_orbit):
        #         self.Nlm = int(block_size / 2)
        #     else:
        #         self.Nlm = block_size
        # self.l = int((self.Nlm - 1) / 2)

    def run_HI(self, calc_off_diag=False, remove_CF=False, zerotemp=False, called_CF_corr=False, lad_bs=False):
        """
        Runs the Hubbard-I approximation and returns the local Green's functions G_at and G_Gamma.

        Parameters:
        -----------
        calc_off_diag : bool, optional
            Whether to calculate off-diagonal elements. Defaults to False.
        remove_CF : bool, optional
            Whether to remove crystal fields. Defaults to False.
        zerotemp : bool, optional
            If set to True, calculations are done at zero temperature. Otherwise, finite temperature is used. Defaults
            to False.
        called_CF_corr : bool, optional
            Indicates if crystal field corrections were applied. Defaults to False. (This parameter seems unused in
            the method.)
        lad_bs : bool, optional
            If set to True, it assumes ladder-based solver. Defaults to False.

        Returns:
        --------
        G_at : numpy.array
            Local Green's function for an atom.
        G_Gamma : list of numpy.array
            List containing Green's functions for different energy levels.

        Notes:
        ------
        This method calculates the atomic Green's function using the Hubbard-I approximation. The method
        assumes full U interaction matrix and supports spin-orbit considerations. It calls an external solver
        `gf_hi_fullu_int` to compute the main results and then formats them appropriately.
        """

        nlms = self.Nlm * 2  # only one block in this case!
        nlm = self.Nlm
        ummss = np.zeros((nlms, nlms))

        if self.U_int:
            ur, umn, ujmn = self.__set_umatrix(U=self.U_int, J=self.J_hund, T=self.T)
            U = np.zeros((nlms, nlms, nlms, nlms))
            ummss[0:nlm, 0:nlm] = ujmn
            ummss[nlm:nlms, nlm:nlms] = ujmn
            ummss[0:nlm, nlm:nlms] = umn
            ummss[nlm:nlms, 0:nlm] = umn
            for is1 in range(2):
                for is2 in range(2):
                    U[is1 * nlm:(is1 + 1) * nlm, is2 * nlm:(is2 + 1) * nlm, is1 * nlm:(is1 + 1) * nlm,
                    is2 * nlm:(is2 + 1) * nlm] = ur
        else:
            U = self.u4ind
            for i in range(nlms):
                for j in range(nlms):
                    ummss[i, j] = U[i, j, i, j] - U[i, j, j, i]

        # Initialization
        if (self.UseSpinOrbit):
            nlms = 2 * nlm
        else:
            nlms = nlm

        iwmax = int(len(self.G_iw.mesh))
        gf0 = np.zeros((nlm * self.n_spin, nlm * self.n_spin, iwmax), dtype=np.cdouble, order='F')
        tail0 = np.zeros((self.Nmoments, nlm * self.n_spin, nlm * self.n_spin), dtype=np.cdouble, order='F')
        gf = np.zeros((self.n_lev, self.n_lev, nlm * self.n_spin, nlm * self.n_spin, iwmax), dtype=np.cdouble,
                      order='F')
        tail = np.zeros((self.n_lev, self.n_lev, self.Nmoments, nlm * self.n_spin, nlm * self.n_spin),
                        dtype=np.cdouble, order='F')

        N = self.ealmat.shape[0]
        if zerotemp:
            temp = 0.0
        else:
            temp = 1.0 / self.beta
        if lad_bs:
            lad_op = self.lad_op
        else:
            lad_op = np.zeros([N, N], np.complex_)
        if not isinstance(self.st_bas, np.ndarray):
            if self.n_bas == None:
                MultNat = 1
                self.n_bas = self.n_lev  # just initializing st_bas dimensions with something to have a proper array
            else:
                MultNat = factorial(N) / factorial(N - self.gs_occ) / factorial(self.gs_occ)
            self.st_bas = np.zeros([MultNat, self.n_bas], np.complex_)
        else:
            MultNat = self.st_bas.shape[0]
            self.n_bas = self.st_bas.shape[1]

        M = [x for x in self.G_iw.mesh]
        self.zmsb = np.array([x for x in M], np.complex_)
        self.ovlmat = np.zeros((self.n_bas, self.n_lev), dtype=np.complex_, order='F')

        self.__save_eal("eal_in_HubI_interact")

        mpi.report("n_lev = %s" % (self.n_lev))
        mpi.report("n_bas = %s" % (self.n_lev))
        mpi.report("calc_off_diag = %s" % (calc_off_diag))
        mpi.report("remove_CF = %s" % (remove_CF))
        mpi.report("nlms  : %s" % (nlms))
        mpi.report("temp  : %s" % (temp))
        mpi.report("lad_bs  : %s" % (lad_bs))
        self.atocc, self.atmag, Z = gf_hi_fullu_int(gf0=gf0, tail0=tail0,
                                                    gf=gf, tail=tail,
                                                    iwmax=iwmax,
                                                    e0f=self.ealmat, u=U,
                                                    ummss=ummss,
                                                    zmsb=self.zmsb, nlm=self.Nlm,
                                                    nmom=self.Nmoments,
                                                    ns=self.n_spin,
                                                    temp=temp,
                                                    verbosity=self.verbosity,
                                                    n_lev=self.n_lev,
                                                    calc_off_diag=calc_off_diag,
                                                    remove_cf=remove_CF,
                                                    ladbs=lad_bs,
                                                    ladop=lad_op,
                                                    nbas=self.n_bas, mnat=MultNat,
                                                    calcovl=self.CalcOvl,
                                                    stbas=self.st_bas,
                                                    ovlmat=self.ovlmat)

        mpi.report("Atomic occupancy in Hubbard I Solver  : %s" % self.atocc)
        mpi.report("Atomic magn. mom. in Hubbard I Solver : %s" % self.atmag)
        mpi.report("Atomic Z  : %s" % Z)
        mpi.report("Nmoments  : %s" % (self.Nmoments))

        if lad_bs:
            return
        elif self.CalcOvl:
            self.ovlmat = self.ovlmat

        G_at = self.__copy_to_GF(gf0, nlms)
        # G_at.save('Gat')
        G_Gamma = [[] for i in range(self.n_lev)]
        for ilev in range(self.n_lev):
            if not calc_off_diag:
                G_Gamma[ilev].append(self.__copy_to_GF(gf[ilev, ilev], nlms))
            else:
                for ilev1 in range(self.n_lev):
                    G_Gamma[ilev].append(self.__copy_to_GF(gf[ilev, ilev1], nlms))

        return G_at, G_Gamma

    def calc_Sig_lev(self, calc_off_diag=False, remove_CF=False):
        """
        Calculate self-energies for the first n_lev atomic levels.

        Parameters:
        -----------
        calc_off_diag : bool, optional
            Whether to calculate off-diagonal elements. Defaults to False.
        remove_CF : bool, optional
            Whether to remove crystal fields. Defaults to False.

        Returns:
        --------
        Sig_lev : list of lists of Green's functions
            The computed self-energies for each energy level.

        Notes:
        ------
        This method calculates self-energies by considering the difference between the Green's function for individual
        atomic levels and the average atomic Green's function. It supports calculations for both diagonal and
        off-diagonal elements of the Green's function matrix.
        """

        mpi.report("\nStarting calculations of atomic GF for GS multiplet levels")

        G_at, G_Gamma = self.run_HI(calc_off_diag=calc_off_diag, remove_CF=remove_CF)

        # DEBUG print begins
        GMat_pos={}
        GMat_neg={}
        for i0 in range(self.n_lev):
            for i1 in range(self.n_lev):
                mesh=np.array([w.real for w in G_at['ud'].mesh])
                iw0=int(mesh.shape[0]/2)
                iw0_n=int(mesh.shape[0]/2-1)
                GMat_pos['%s_%s'%(i0,i1)]=G_Gamma[i0][i1]['ud'].data[iw0,:,:]
                GMat_neg['%s_%s'%(i0,i1)]=G_Gamma[i0][i1]['ud'].data[iw0_n,:,:]
        if mpi.is_master_node():
            ar=HDFArchive('GMat.h5','a')
            ar['GMat_pos']=GMat_pos
            ar['GMat_neg']=GMat_neg
            del ar        
        # END DEBUG

        Sig_lev = [[] for i in range(self.n_lev)]
        Gat_inv = inverse(G_at)
        Gat_av = G_at.copy()
        G_GamGam = []
        # copy diagonal elements of Gmat to list G_GamGam
        # note: if calc_off_diag=False their are in G_Gamma[ilev][0]
        if not calc_off_diag:
            for ilev in range(self.n_lev): G_GamGam.append(G_Gamma[ilev][0])
        else:
            for ilev in range(self.n_lev): G_GamGam.append(G_Gamma[ilev][ilev])

        for ilev in range(self.n_lev):
            if (ilev == 0):
                Gat_av << G_GamGam[ilev]
            else:
                Gat_av += G_GamGam[ilev]
        Gat_av /= self.n_lev

        for ilev in range(self.n_lev):
            if not calc_off_diag:
                Sig_lev[ilev].append(Gat_inv * (G_GamGam[ilev] - Gat_av) * Gat_inv)
            else:
                for ilev1 in range(self.n_lev):
                    if ilev == ilev1:
                        Sig_lev[ilev].append(Gat_inv * (G_GamGam[ilev] - Gat_av) * Gat_inv)
                    else:
                        Sig_lev[ilev].append(Gat_inv * G_Gamma[ilev1][ilev] * Gat_inv)

        self.G_Gamma = G_Gamma
        return Sig_lev

    def __copy_to_GF(self, gf, nlmtot):
        """ Put GF form hub-I solver into GF class"""
        # transfer the data to the GF class:
        M = {}
        isp = -1
        for a, al in self.gf_struct:
            isp += 1
            M[a] = np.array(gf[isp * nlmtot:(isp + 1) * nlmtot, isp * nlmtot:(isp + 1) * nlmtot, :]).transpose(2, 0,
                                                                                                               1).copy()
        glist = []
        name_list = [block for block, block_size in self.gf_struct]
        for block, block_size in self.gf_struct:
            glist.append(
                GfImFreq(beta=self.beta, n_points=self.n_iomega, target_shape=(len(block_size), len(block_size))))

        G = BlockGf(name_list=name_list, block_list=glist)
        self.__copy_Gf(G, M)

        # glist = lambda: [GfImFreq(indices=al, beta=self.beta, n_points=self.n_omega) for a, al in self.gf_struct]
        # G = BlockGf(name_list=self.a_list, block_list=glist(), make_copies=False)
        # self.__copy_Gf(G, M)

        return G

    def __set_umatrix(self, U, J, T=None):
        if self.l > 0:
            Umat = U_matrix(l=self.l, U_int=U, J_hund=J, T=T)
            U, Up = reduce_4index_to_2index(Umat)
        else:
            Umat = np.zeros((1, 1, 1, 1), dtype=float)
            Umat[0, 0, 0, 0] = U
            U = np.zeros((1, 1), dtype=float)
            Up = np.zeros((1, 1), dtype=float)
            Up[0, 0] = U

        return Umat, Up, U

    def set_ud_levels(self, eal, rmat=None,rmat_time_inv=0):
        """
        Set atomic levels into a single `ud` block.

        This method is used to combine the spin-up (`up`) and spin-down (`down`) atomic levels
        into a single `ud` block representing both spin channels. If the input already contains
        the 'ud' key, atomic levels are directly set using this information. Otherwise, atomic levels
        are constructed by stacking the 'up' and 'down' matrices.

        Parameters:
        -----------
        eal : dict
            Dictionary containing atomic levels. The keys can be 'ud', 'up', and 'down',
            representing spin-up-down combined, spin-up, and spin-down atomic levels, respectively.
            Each key maps to a matrix representing atomic energy levels for that spin.
        rmat : numpy ndarray
             rotation matrix in spin-orbital basis
        rmat_time_inv: integer
             =1 if the time reversal operation is to be applied together with rotation

        Notes:
        ------
        The constructed `ud` block will have twice the rows and columns of the 'up' or 'down' block.
        The top-left and bottom-right quarters of the `ud` matrix correspond to 'up' and 'down' atomic
        levels, respectively.
        Rotation is applied to 'ud' block if the rotation matrix rmat is submitted

        """

        ealmat = {}
        if 'ud' in eal:
            ealmat['ud']=eal['ud'].copy()
        else:
            Nlm = eal['up'].shape[0]
            ealmat['ud'] = np.zeros((Nlm * 2, Nlm * 2), np.complex)
            ealmat['ud'][0:Nlm, 0:Nlm] = eal['up']
            ealmat['ud'][Nlm:2 * Nlm, Nlm:2 * Nlm] = eal['down']
        nlms=ealmat['ud'].shape[0]
        nlm=int(nlms/2)
        if rmat.shape[0]==nlm:
            # only orbital-space rotation is given (no SO)
            rmat_tmp=rmat.copy()
            rmat=np.zeros((nlms,nlms),complex)
            rmat[0:nlm,0:nlm]=rmat_tmp
            rmat[nlm:nlms,nlm:nlms]=rmat_tmp
        if isinstance(rmat,np.ndarray):
            if rmat_time_inv==1:
                ealmat['ud']=np.dot(rmat.conjugate(),np.dot(ealmat['ud'].transpose(),rmat.transpose()))
            else:
                ealmat['ud']=np.dot(rmat,np.dot(ealmat['ud'],rmat.transpose().conjugate()))
        self.set_atomic_levels(eal=ealmat)


    def __save_eal(self,Filename):
        if mpi.is_master_node():
            f=open(Filename,'a')
            f.write('\neff. atomic levels\n')
            N=self.ealmat.shape[0]
            for i in range(N):
                for j in range(N):
                    f.write("%10.6f %10.6f   "%(self.ealmat[i,j].real,self.ealmat[i,j].imag))
                f.write("\n")
            f.close()


    def __copy_Gf(self, G, data):
        """ Copies data and tail to Gf object GF """
        for s, g in G:
            g.data[:, :, :] = data[s][:, :, :]
