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
from numpy import *
import numpy as np
from triqs.gf import *
import triqs.utility.mpi as mpi
from triqs_dft_tools.symmetry import *
from h5 import *
from MagInt.HubbardI_interact import *
from MagInt.utils import *
from MagInt.Shells import *
from MagInt.Read_kp_sym import *
from hubbardI.hubbard_I import fock_space_rmat
import sys


class MagInteract:
    """ Implements the force-theorem Hubbard-I (FT-HI) approach to intersite
    exchange interactions in correlated insulators.
    The formalism is given in L. V. Pourovskii Phys. Rev. B 94, 115117 (2016)
    This Python-3 version is based on TRIQS library """

    def __init__(self, general_par, magint_par, label_corrsite, GSM=None, SK=None, S_HI=None, GF_tmat=None):
        """
        Constructor for initializing magnetic interaction calculations.

        Parameters:
        -----------
        general_par : dict
            Dictionary containing general parameters, e.g., 'dft_exec' specifying the DFT execution type and 'filename'.
        magint_par : dict
            Dictionary containing parameters specific to magnetic interactions, e.g., 'atom1', 'atom2', and 'n_shells'.
        label_corrsite : dict
            Label for the correlated sites.
        GSM : object, optional
            Green's function mesh. None by default.
        SK : object, optional
            Object containing structural information. None by default.
        S_HI : object, optional
            Object containing information on the Hamiltonian interaction structure. None by default.
        GF_tmat : list or None, optional
            List of transformation matrices to apply to intersite GF at each k-point, if required. None by default.
        """

        self.dft_exec = general_par['dft_exec']
        self.SK = SK
        self.lab_cr = label_corrsite
        self.Interactions = ((magint_par['atom1'], magint_par['atom2']),)
        self.max_n_shells = magint_par['n_shells']
        self.fname = general_par['filename']
        self.GSM = GSM
        self.S_HI = S_HI
        self.interact_types = {}  # set()
        self.interact_sites = {}
        self.interact_basis = {}
        self.calc_all_sites = magint_par['calc_all_sites']

        for intr in self.Interactions:
            for type in intr:
                # self.interact_types.add(type)
                self.interact_types[type] = None
                if type not in self.interact_sites:
                    self.interact_sites[type] = []
                    self.interact_basis[type] = 0
                    count = 0
                    for icrsh in SK.corr_to_inequiv:
                        # if SK.corr_to_inequiv[icrsh] == self.lab_cr[type]:
                        if icrsh == self.lab_cr[type]:
                            # self.interact_sites[type].append(icrsh)
                            self.interact_sites[type].append((icrsh, count))
                            self.interact_basis[type] += 1
                            count += 1

        if GSM and SK and S_HI:
            self.Mesh = self.SK.Sigma_imp_iw[0].mesh
            mpi.report("-" * 40)
            mpi.report("\nAtomic types between which interactions will be calculated:")
            mpi.report("\nType   GSM degeneracy    Ineq. Corr. Site   Equiv. Corr. Sites        Num. Basis Vec.")
            for type in self.interact_types:
                mpi.report("%s      %s             %s                  %s               %s" % (
                    type, self.GSM[type], self.lab_cr[type], self.interact_sites[type], self.interact_basis[type]))
        elif SK:
            mpi.report("\n GSM is not set, MagInt calculations impossible. Only effective hopping calcs. available")
        else:
            "SK needs to be set!"
            assert 0
        # DEBUG: GF_tmat is a list of transformation matrices t
        if GF_tmat:
            # transformation matrices to apply to intersite GF at each k-point
            mpi.report("\n Inter-site GF will be transformed using provided GF_tmat")
            self.trans_GF = True
            self.GF_tmat = GF_tmat
        else:
            self.trans_GF = False
        # end DEBUG

    def rot_stbas(self, StBas, n_el, nlm, rot_mat_in, rot_time_rev):
        """
        Rotate standard basis 
 
        Parameters:
        -----------
        StBas : np.ndarray
            initial standard basis
        n_el : integer
            shell occupancy
        nlm  : integer
            shell orbital degeneracy
        rot_mat_in : np.ndarray
            one-electron rotation matrix 
            of dimension ( nlm*2 x nlm*2 ) if SO included or ( nlm x nlm ) otherwise
        rot_time_rev : integer
            =1 if time reversal is to be applied together with rotation
       
        Returns:
        --------
        StBas_out : np.ndarray
            transformed Standard Basis
        """
        nlms = 2 * nlm
        if rot_mat_in.shape[0] == nlm:
            # rotmat without SO, extending it to both spin blocks
            rot_mat=np.zeros((nlms,nlms),complex)
            rot_mat[0:nlm,0:nlm]=rot_mat_in
            rot_mat[nlm:nlms,nlm:nlms]=rot_mat_in
        else:
            rot_mat=rot_mat_in
        vlen = int(factorial(nlms) / factorial(n_el) / factorial(nlms - n_el))
        fs_rmat = np.zeros((vlen, vlen), dtype=np.cdouble, order='F')
        if rot_time_rev == 1:
            # time reversal included
            fock_space_rmat(fs_rmat=fs_rmat, rm=rot_mat.transpose(), nat=n_el, nso=nlms, n=vlen)
        else:
            fock_space_rmat(fs_rmat=fs_rmat, rm=rot_mat, nat=n_el, nso=nlms, n=vlen)
        N_lev = StBas.shape[1]
        StBas_out = np.zeros((vlen, N_lev), complex)
        for i in range(N_lev):
            if rot_time_rev == 1:
                # time reversal included
                StBas_out[:, i] = np.dot(fs_rmat, StBas[:, i].conjugate())
            else:
                StBas_out[:, i] = np.dot(fs_rmat, StBas[:, i])
        return StBas_out

    def calc_MagInt(self, general_par, S_INT, one_per_shell=False,
                    calc_off_diag=False, remove_CF=False, R_arr=None, R_arr_site=None, Sigma_fix=None):
        """
        Perform calculation of magnetic pair interactions by considering the
        intersite Green's function and the Hubbard-I self-energy.

        The function calculates magnetic pair interactions by first reading the real-space
        intersite vectors R, generating k-mesh in the full Brillouin zone, and then calculating
        the intersite Green's function G_R. The interactions are then computed with the Hubbard-I
        self-energy. The interactions are stored and printed in meV.

        Parameters:
        -----------
        general_par : dict
            General parameters required for the calculation.
        S_INT : dict
            Hubbard-I self-energy.
        one_per_shell : bool, optional
            Calculate only one interaction for each coordinate shell. Default is False.
        calc_off_diag : bool, optional
            Calculate off-diagonal interactions. Default is False.
        remove_CF : bool, optional
            Remove crystal field. Default is False.
        R_arr : np.ndarray or None, optional
            Real-space lattice vectors connecting correlated sites. Default is None.
        R_arr_site : np.ndarray or None, optional
            Sites corresponding to the real-space vectors in `R_arr`. Default is None.
        Sigma_fix : np.ndarray or None, optional
            Fixed self-energy input. If provided, Sigma will be set to this value. Default is None.

        Returns:
        --------
        V : dict
            Calculated magnetic interactions.
        V_4ind : dict
            Full 4-index magnetic interactions if `calc_off_diag` is True.
        int_pairs : dict
            Interaction pairs for the calculated interactions.

        Notes:
        ------
        For Wien2k executions, the lattice vectors are read from the ".outputnn" file.
        Local rotations are suppressed during the calculation.
        """

        self.int_pairs = {}

        if general_par['dft_exec'] == 'Wien2k':
            self.latt = zeros([3])
            f = open(self.fname + ".outputnn")
            head = [f.readline() for i in range(4)]
            width = 10
            #rela_line = head[3].strip()
            rela_line = head[3]
            mpi.report("%s\n"%rela_line)
            self.latt = np.array([float(rela_line[j:j+width].strip()) for j in range(0, 3*width, width)])
            if general_par['verbosity'] > 0:
                mpi.report("\nA = %10.6f    B = %10.6f    C = %10.6f" % (self.latt[0], self.latt[1], self.latt[2]))
            f.close()

        for intr in self.Interactions:
            if general_par['verbosity'] > 0:
                mpi.report("\n----INTERACTION TYPE: %s-%s" % (intr[0], intr[1]))

            self.int_pairs[intr] = shells(general_par, self.fname, intr, self.interact_sites, self.interact_basis,
                                          max_n_shells=self.max_n_shells, print_vecs='Basis', R_arr=R_arr,
                                          R_arr_site=R_arr_site, calc_all_sites=self.calc_all_sites)

            self.int_pairs[intr].calc_d_origin()

        # Generate the k-mesh
        self.KP = read_kp_sym(general_par)

        # Calculate delta_Sigma for each interacting type and site in the global frame
        Sig_lev = {}
        for type in self.interact_types:
            Sig_lev[type] = {}
            # Sig_tmp = S_INT[type].calc_Sig_lev(calc_off_diag, remove_CF)
            # if Sigma_fix:
            #    mpi.report("\nSigma is set to input Sigma_fix\n")
            #    for S in Sig_tmp:
            #        for sig in S:
            #            sig << Sigma_fix
            for val in self.interact_sites[type]:
                ieq = val[1]
                Sig_lev[type][ieq] = S_INT[type][ieq].calc_Sig_lev(calc_off_diag, remove_CF)

        # DEBUG: print Sig_lev
        # for type in self.interact_types:
        #    for val in self.interact_sites[type]:
        #        ieq = val[1]
        #        print('Row number: ',len(Sig_lev[type][ieq]))
        #        for i in range(len(Sig_lev[type][ieq])): print(len(Sig_lev[type][ieq][i]))
        #        Sigma=Sig_lev[type][ieq][0][0]['ud']
        #        mesh=np.array([w.real for w in Sigma.mesh])
        #        iw0=int(mesh.shape[0]/2)
        #        mat=Sigma.data[iw0,:,:]
        #        if mpi.is_master_node(): print_arr(mat,log='\nSig_lev site=%s'%ieq)
        #        if ieq==0:
        #            rmat=self.SK.rot_mat
        #            mat1=np.dot(rmat[1],np.dot(mat,rmat[1].transpose().conjugate()))
        #            if mpi.is_master_node(): print_arr(mat1,log='\nSig_lev after rotation')

        # assert 0
        # END DEBUG

        # Collect interactions to be calculated into dicts
        self.sig_store = Sig_lev
        G_off = {}
        G_off_inv = {}
        for intr in self.Interactions:
            pair = self.int_pairs[intr]
            for key0 in pair.R_vec:
                type0 = key0[0][0]
                type1 = key0[1][0]
                mul0 = key0[0][1]
                mul1 = key0[1][1]
                for key1 in pair.R_vec[key0]:

                    key = "%s_%s_%s_%s_%s" % (type0, mul0, type1, mul1, key1)

                    # DEBUG
                    name_list = [a for a, al in S_INT[type0][mul0].gf_struct]
                    g_iw_list = []
                    for a, al in S_INT[type0][mul0].gf_struct:
                        g_iw_list.append(GfImFreq(mesh=self.Mesh, target_shape=(len(al), len(al))))
                    G_off[key] = BlockGf(name_list=name_list, block_list=g_iw_list)
                    G_off_inv[key] = G_off[key].copy()

                    if one_per_shell:
                        # calculate only one interaction for each coor. shell
                        break

        # calculate inter-site GF
        self.inter_site_GF(general_par, G_off, G_off_inv)

        # calculate the interactions
        mpi.report('\n\n ---------- CALCULATING MAGNETIC INTERACTIONS:\n')
        V = {}
        V_4ind = {}
        for intr in self.Interactions:
            pair = self.int_pairs[intr]
            for key0 in pair.R_vec:
                (type0, mul0), (type1, mul1) = key0
                for key1 in pair.R_vec[key0]:

                    # site0 = self.interact_sites[type0][0]
                    # site1 = self.interact_sites[type1][0]
                    dim0 = self.GSM[type0]
                    dim1 = self.GSM[type1]
                    Sig0 = Sig_lev[type0]
                    Sig1 = Sig_lev[type1]

                    key = "%s_%s_%s_%s_%s" % (type0, mul0, type1, mul1, key1)
                    mpi.report("%s" % (key))

                    # print("\nDEBUGGONE")
                    # print(G_off[key]['ud'][0, 0])
                    # print(Sig0[mul0][0][0]['ud'][0, 0])
                    # print(Sig1[mul1][0][0]['ud'][0, 0])

                    if (not key in G_off): continue
                    if calc_off_diag:
                        V_4ind[key] = zeros([dim0, dim0, dim1, dim1], complex_)
                        for i in range(dim0):
                            for j in range(dim0):
                                for k in range(dim1):
                                    for l in range(dim1):
                                        # dens = (G_off[key] * Sig1[site1][l][k] * G_off_inv[key] * Sig0[site0][i][j]).density()
                                        dens = (G_off[key] * Sig1[mul1][l][k] * G_off_inv[key] * Sig0[mul0][i][
                                            j]).density()
                                        for s, g in G_off[key]:
                                            V_4ind[key][i, j, l, k] += np.trace(dens[s])

                        V[key] = zeros([dim0, dim1], complex_)
                        for i in range(dim0):
                            for j in range(dim1):
                                V[key][i, j] = V_4ind[key][i, i, j, j]
                    else:
                        V[key] = zeros([dim0, dim1], complex_)
                        for i in range(dim0):
                            for j in range(dim1):
                                # V[key][i, j] = (G_off[key] * Sig0[site0][i][0] * G_off_inv[key] * Sig1[site1][j][
                                V[key][i, j] = (G_off[key] * Sig0[mul0][i][0] * G_off_inv[key] * Sig1[mul1][j][
                                    0]).total_density()

        mpi.report("\n\n ---------- MAGNETIC INTERACTIONS in meV")
        for intr in self.Interactions:
            pair = self.int_pairs[intr]
            for key0 in pair.R_vec:
                (type0, mul0), (type1, mul1) = key0
                dim0 = self.GSM[type0]
                dim1 = self.GSM[type1]
                mpi.report("\n\n ----- INTERACTION: %s - %s" % (type0, type1))
                for i in range(pair.nshells):
                    mpi.report("\nShell = %s" % (i))
                    mpi.report("Coor_num = %s" % (pair.coor_num[key0][i]))
                    for i1 in range(pair.coor_num[key0][i]):
                        key = "%s_%s_%s_%s_%s_%s" % (type0, mul0, type1, mul1, i, i1)
                        if (key in V):
                            mpi.report("\n key= %s\n" % (key))
                            R_key = "%s_%s" % (i, i1)
                            mpi.report(
                                "\n R = %s, sites: 0 - %s\n" % (pair.R_vec[key0][R_key], pair.R_site[key0][R_key]))
                            for k in range(dim0):
                                st = ''
                                for l in range(dim1):
                                    st += "  %8.3f" % (V[key][k, l].real * 1000.0)
                                mpi.report(st)

        return V, V_4ind, self.int_pairs

    def inter_site_GF(self, general_par, G_off, G_off_inv):
        """
        Calculates inter-site Green's Functions (GF) for all required pair vectors.

        The function computes the inter-site GF based on the Bloch's theorem and updates the `G_off` dictionary with
        the Fourier-transformed GF values. This function also handles the inversion symmetry to update the `G_off_inv`
        dictionary.

        Parameters:
        -----------
        general_par : dict
            Dictionary of general parameters. Must include 'dft_exec' which specifies the DFT execution type
            (either 'Wien2k' or 'Vasp').
        G_off : dict
            Dictionary to store the Fourier-transformed inter-site GF. The keys are constructed from type0, type1,
            and R_key.
        G_off_inv : dict
            Dictionary to store the Fourier-transformed inter-site GF considering inversion symmetry. It has the
            same structure as G_off.

        Returns:
        --------
        None.
            The method updates the dictionaries `G_off` and `G_off_inv` in-place.

        Raises:
        -------
        AssertionError
            If the given block names in G_off are neither 'up/down' nor 'ud'.
            If there's an inconsistency between direct and inv transformed sites.
            If there's a key error for constructed keys in G_off.
        """

        SK = self.SK
        KP = self.KP

        # temporary storage
        GF_tmp = {}
        for int in self.Interactions:
            type0 = int[0]
            name_list = [a for a, al in self.S_HI[type0].gf_struct]
            g_iw_list = []
            for a, al in self.S_HI[type0].gf_struct:
                g_iw_list.append(GfImFreq(mesh=self.Mesh, target_shape=(len(al), len(al))))
            GF_tmp[type0] = BlockGf(name_list=name_list, block_list=g_iw_list)

        for ik in mpi.slice_array(np.array(range(SK.n_k))):
            mpi.report("  IK = %s" % (ik))

            GF_k_lat = SK.lattice_gf(ik)

            for int in self.Interactions:

                int_pair = self.int_pairs[int]
                type0 = int[0]
                # type1 = int[1]

                G_t = GF_tmp[type0]
                G_t_inv = GF_tmp[type0].copy()

                # fsite0 = self.interact_sites[type0][0]
                # fsite1 = self.interact_sites[type1][0]

                if ('up' in G_t.indices) and ('down' in G_t.indices):
                    Nlm = G_t['up'].data.shape[1]
                    SOblocks = False
                elif 'ud' in G_t.indices:
                    SOblocks = True
                else:
                    assert 0, 'Block names in G_off should be up/down or ud!'

                for key0 in int_pair.T_vec:
                    (type0, mul0), (type1, mul1) = key0
                    ish0 = self.interact_sites[type0][0][0]
                    ish1 = self.interact_sites[type1][0][0]


                    d_origins_init = int_pair.d_origin[key0]

                    for sig, gf in G_t:
                        isp = SK.spin_names_to_ind[SK.SO][sig]
                        n_orb = SK.n_orbitals[ik, isp]
                        projmat0 = SK.proj_mat[ik][isp][ish0][:, 0:n_orb]
                        projmat1 = SK.proj_mat[ik][isp][ish1][:, 0:n_orb]
                        gf.from_L_G_R(projmat0, GF_k_lat[sig], projmat1.conjugate().transpose())

                    for sig, gf in G_t_inv:
                        isp = SK.spin_names_to_ind[SK.SO][sig]
                        n_orb = SK.n_orbitals[ik, isp]
                        projmat0 = SK.proj_mat[ik][isp][ish0][:, 0:n_orb]
                        projmat1 = SK.proj_mat[ik][isp][ish1][:, 0:n_orb]
                        gf.from_L_G_R(projmat1, GF_k_lat[sig], projmat0.conjugate().transpose())

                    if self.trans_GF:
                        for sig, gf in G_t:
                            gf.from_L_G_R(self.GF_tmat[ik], gf, self.GF_tmat[ik])
                        for sig, gf in G_t_inv:
                            gf.from_L_G_R(self.GF_tmat[ik], gf, self.GF_tmat[ik])

                    if general_par['dft_exec'] == 'Wien2k':
                        array_kpt = len(KP.k_fbz[ik])
                    elif general_par['dft_exec'] == 'Vasp':
                        array_kpt = 1

                    for ikp in range(array_kpt):

                        if general_par['dft_exec'] == 'Wien2k':
                            syop = KP.syop_fbz[ik][ikp]
                            K_shift = KP.K_shift[ik][ikp]

                            # TODO DEBUG!
                            # if (SK.symmcorr.orb_map[syop][ish0] != fsite0): continue
                            if general_par['use_symmetries']:

                                M, ish0_n, ish1_n = self.__rot_GF_(G_t, syop, ish0, ish1, K_shift, d_origins_init)
                                M_inv, dumm1_n, dumm0_n = self.__rot_GF_(G_t_inv, syop, ish1, ish0, K_shift,
                                                                         -d_origins_init)
                                if (ish0_n != dumm0_n) or (ish1_n != dumm1_n):
                                    print('Inconsistency between direct/inv transformed sites')
                                    print('for ik=%s ikp=%s, symop=%' % (ik, ikp, syop))
                                    assert 0

                            else:
                                M = G_t;
                                M_inv = G_t_inv

                            kvec = KP.k_fbz[ik][ikp]

                        elif general_par['dft_exec'] == 'Vasp':

                            kvec = KP.k_fbz[0][ik]
                            M = G_t
                            M_inv = G_t_inv

                        for R_key, R in int_pair.T_vec[key0].items():

                            key = "%s_%s_%s_%s_%s" % (type0, mul0, type1, mul1, R_key)
                            if not key in G_off:
                                print('Key error for key ', key)
                                assert 0

                            # TODO VASP Check for supercell
                            # The phase parameter is set to 1 as per default and is meant for debugging
                            phase = - 2.0 * pi * dot(kvec, R) * general_par['phase']
                            exp_fac = exp(complex(0.0, phase))

                            if SOblocks:
                                # copy into sig=ud block of  G_off
                                for sig, gf in M:
                                    G_off[key][sig] += M[sig] * exp_fac * KP.weight
                                    G_off_inv[key][sig] += M_inv[sig] * np.conjugate(exp_fac) * KP.weight
                            else:
                                # collect up/down contirbuitons
                                G_off[key]['ud'][0:Nlm, 0:Nlm] += M['up'] * exp_fac * KP.weight
                                G_off[key]['ud'][Nlm:2 * Nlm, Nlm:2 * Nlm] += M['down'] * exp_fac * KP.weight
                                G_off_inv[key]['ud'][0:Nlm, 0:Nlm] += M_inv['up'] * np.conjugate(
                                    exp_fac) * KP.weight
                                G_off_inv[key]['ud'][Nlm:2 * Nlm, Nlm:2 * Nlm] += M_inv['down'] * np.conjugate(
                                    exp_fac) * KP.weight

        for key in G_off:
            G_off[key] = mpi.all_reduce(mpi.world, G_off[key], lambda x, y: x + y)
            G_off_inv[key] = mpi.all_reduce(mpi.world, G_off_inv[key], lambda x, y: x + y)

    def __rot_GF_(self, GF_k, syop, iorb0, iorb1, K_shift, dtau):
        """
        Apply a symmetry operation to the Hamiltonian in k-space.

        Parameters:
        -----------
        self : object
            The instance of the class.
        H_k : ndarray
            The Hamiltonian matrix in k-space.
        syop : int
            The index of the symmetry operation to apply.
        iorb0 : int
            The index of the initial orbital.
        iorb1 : int
            The index of the final orbital.
        K_shift : ndarray, optional
            A vector representing the shift in k-space due to the symmetry operation. Default is None.
        dtau : ndarray, optional
            A vector representing the displacement in real space. Default is None.

        """
        Symm = self.SK.symmcorr
        GF_rot = GF_k.copy()
        GF_cp = GF_k.copy()
        # print '\n',syop,iorb0,iorb1
        # print dtau
        iorb_n0 = Symm.orb_map[syop][iorb0]
        iorb_n1 = Symm.orb_map[syop][iorb1]
        phase = dot(K_shift, dtau)
        if (Symm.time_inv[syop]): GF_cp << GF_k.transpose()
        for sig, gf in GF_rot:
            # print(Symm.mat[syop][iorb0])
            gf.from_L_G_R(Symm.mat[syop][iorb0], GF_cp[sig], Symm.mat[syop][iorb1].conjugate().transpose())
            # gf *= exp(complex(0.0, -2.0 * pi * phase))
        return GF_rot, iorb_n0, iorb_n1

    def __rotsig_(self, GF, icrsh, fsite):
        """ rotates sigma (or GF) into equivalent site icrsh"""
        Symm = self.SK.symmcorr
        for syop in range(Symm.n_symm):
            if (Symm.map[syop][fsite] == icrsh):
                if (Symm.time_inv[syop]): GF <<= GF.transpose()
                for sig, gf in GF: gf.from_L_G_R(Symm.mat[syop][fsite], gf,
                                                 Symm.mat[syop][fsite].conjugate().transpose())
                return

    def __rot_H_(self, H_k, syop, iorb0, iorb1, K_shift=None, dtau=None):
        """apply symmetry operation syop to H_k"""
        Symm = self.SK.symmcorr
        iorb_n0 = Symm.orb_map[syop][iorb0]
        iorb_n1 = Symm.orb_map[syop][iorb1]
        if K_shift and dtau: phase = dot(K_shift, dtau)
        if Symm.time_inv[syop] == 0:
            H_new = np.dot(Symm.mat[syop][iorb0], np.dot(H_k, Symm.mat[syop][iorb1].conjugate().transpose()))
        else:
            H_new = np.dot(Symm.mat[syop][iorb0],
                           np.dot(H_k.conjugate(), Symm.mat[syop][iorb1].conjugate().transpose()))
        return H_new, iorb_n0, iorb_n1

    def __print_gf_(self, GF, log='\nGF:'):
        """
        Print the Green's Function at the first Matsubara frequency.

        Parameters:
        -----------
        self : object
            The instance of the class.
        GF : BlockGf
            The Green's Function to be printed.
        log : str, optional
            A log message to display before printing the Green's Function. Default is '\nGF:'.

        This function prints the real and imaginary parts of the Green's Function at the first Matsubara frequency.

        """

        mpi.report("%s\n\n\n" % log)
        for sig, gf in GF:
            dd = gf.data.shape
            nw = dd[0]
            r_part = []
            i_part = []
            for i in range(dd[1]):
                str_r = ''
                str_i = ''
                for k in range(dd[2]):
                    str_r += " %10.7f" % (gf.data[nw / 2, i, k].real)
                    str_i += " %10.7f" % (gf.data[nw / 2, i, k].imag)
                r_part.append(str_r)
                i_part.append(str_i)
        mpi.report("\n----REAL:")
        for r in r_part:  mpi.report(r)
        mpi.report("\n----IMAG:")
        for i in i_part:  mpi.report(i)
        tr = np.trace(gf.data[nw / 2, :, :])
        mpi.report("\n Trace = %12.7f %12.7f" % (tr.real, tr.imag))

    def __print_arr_(self, arr, log='arr'):
        """
        Print a NumPy array with a given label.

        Parameters:
        -----------
        self : object
            The instance of the class.
        arr : numpy.ndarray
            The NumPy array to be printed.
        log : str, optional
            A label to display before printing the array. Default is 'arr'.

        This function prints the elements of a NumPy array with a specified label. If the array contains complex numbers
        , both the real and imaginary parts are displayed.

        """

        m, n = arr.shape
        cmplx = False
        if np.sum(abs(arr.imag)) > 1e-7:
            cmplx = True
        mpi.report("%s" % log)
        for i in range(m):
            st = ''
            for j in range(n):
                if cmplx:
                    st += ' %9.5f %9.5f' % (arr[i, j].real, arr[i, j].imag)
                else:
                    st += ' %9.5f' % arr[i, j].real
            mpi.report("%s" % st)

    def reorder_V(self, V, V4, orl):
        """
        Reorders V4 in accordance with the permutation matrix order.

        Parameters:
        -----------
        self : object
            The instance of the class.
        V : dict
            A dictionary containing complex-valued matrices.
        V4 : dict
            A dictionary containing complex-valued rank-4 tensors.
        orl : list
            A list containing the permutation order.

        This function reorders the rank-4 tensors in `V4` according to the permutation order `orl`.
        The result is stored back in the `V4` dictionary.

        """
        for key, Vmat in V4.items():
            Vc = np.copy(Vmat)
            dim = Vmat.shape
            for j in range(dim[0]):
                for k in range(dim[1]):
                    for l in range(dim[2]):
                        for m in range(dim[3]):
                            Vmat[j, k, l, m] = Vc[orl[j], orl[k], orl[l], orl[m]]
            for j in range(dim[0]):
                for k in range(dim[2]):
                    V[key][j, k] = Vmat[j, j, k, k]

    def trans_V(self, V, V4, S_INT, tol=1e-7, site_list=None):
        """
        Transform V and V4 to a new basis given by transmat.

        Parameters:
        -----------
        self : object
            The instance of the class.
        V : dict
            A dictionary of matrices to be transformed.
        V4 : dict
            A dictionary of 4-index matrices to be transformed.
        transmat0 : np.ndarray or list
            Transformation matrix for the first site. If a list, transmat0[0] is applied to the first site.
        transmat1 : np.ndarray
            Transformation matrix for the second site.
        tol : float, optional
            Tolerance for considering small values as zero. Default is 1e-7.
        site_list : list, optional
            List of two sites for which the transformation is applied. Default is None.

        Returns:
        --------
        None

        Transforms the matrices in V and V4 to a new basis using the provided transformation matrices.

        - If transmat0 is an np.ndarray, it's applied to both sites.
        - If transmat0 is a list, transmat0[0] is applied to the first site, and transmat1 is applied to the second site.

        The transformed matrices are stored back in V and V4.

        """
        Vc = {}
        # TODO Check this!
        # if isinstance(transmat0, type(np.zeros([1, 1]))):
        #     trans0 = transmat0
        #     trans1 = transmat1
        # else:
        #     trans0 = transmat0[0]
        #     trans1 = transmat1[1]

        # DEBUG
        for int in self.Interactions:
            pair = self.int_pairs[int]
            for key0 in pair.R_vec:
                (type0, mul0), (type1, mul1) = key0

                trans0 = S_INT[type0][mul0].ovlmat
                trans1 = S_INT[type1][mul1].ovlmat

                dim0 = trans0.shape[0]
                dim1 = trans0.shape[1]
                dim2 = trans1.shape[0]
                dim3 = trans1.shape[1]
                el = []
                for j in range(dim0):
                    for jj in range(dim1):
                        if abs(trans0[j, jj]) < tol: continue
                        for k in range(dim0):
                            for kk in range(dim1):
                                if abs(trans0[j, jj]) * abs(trans0[k, kk]) < tol: continue
                                for l in range(dim2):
                                    for ll in range(dim3):
                                        if abs(trans0[j, jj]) * abs(trans0[k, kk]) * abs(trans1[l, ll]) < tol: continue
                                        for m in range(dim2):
                                            for mm in range(dim3):
                                                if abs(trans0[j, jj]) * abs(trans0[k, kk]) * abs(trans1[l, ll]) * abs(
                                                        trans1[m, mm]) < tol: continue
                                                val = trans0[j, jj].conjugate() * trans1[l, ll].conjugate() * trans0[
                                                    k, kk] * \
                                                      trans1[m, mm]
                                                el.append([j, k, l, m, jj, kk, ll, mm, val])

                for R_key, R in pair.T_vec[key0].items():
                    # for key, Vmat in V4.items():
                    key = "%s_%s_%s_%s_%s" % (type0, mul0, type1, mul1, R_key)
                    if site_list:
                        if (site_list[0] not in key) or (site_list[1] not in key):
                            continue
                        elif key.index(site_list[0]) > key.index(site_list[1]):
                            continue
                    Vc[key] = np.copy(V4[key])
                    if Vc[key].shape[0] != dim1 or Vc[key].shape[2] != dim3:
                        print('Wrong dimensions of transformation matrix in transmat!')
                        assert 0
                    V4[key] = zeros([dim0, dim0, dim2, dim2], complex_)
                    for itm in el:
                        V4[key][itm[0], itm[1], itm[2], itm[3]] += Vc[key][itm[4], itm[5], itm[6], itm[7]] * itm[8]

                    # for j in range(dim0):
                    #    for k in range(dim0):
                    #        for l in range(dim2):
                    #            for m in range(dim2):
                    #                for jj in range(dim1):
                    #                    for kk in range(dim1):
                    #                        for ll in range(dim3):
                    #                            for mm in range(dim3):
                    #                                #Vmat[j,k,l,m]+=trans0[j,jj]*trans1[l,ll]*trans0inv[kk,k]
                    #                                *trans1inv[mm,m]*Vc[jj,kk,ll,mm]
                    #                                Vmat[j,k,l,m]+=trans0[j,jj].conjugate()*trans1[l,ll].conjugate()
                    #                                *trans0[k,kk]*trans1[m,mm]*Vc[jj,kk,ll,mm]
                    # V4[key] = Vmat
                    V[key] = zeros([dim0, dim2], complex_)
                    for j in range(dim0):
                        for k in range(dim2):
                            V[key][j, k] = V4[key][j, j, k, k]

    def Retrieve_StBas(self, item, h5name='Standard_Basis.h5'):
        """
        Retrieves the standard basis from Standard_Basis.h5.

        Parameters:
        -----------
        self : object
            The instance of the class.
        item : list
            A list containing the group and item names to retrieve from the HDF5 file, e.g., [h5_group, h5_item].
        h5name : str, optional
            The name of the HDF5 file where the standard basis is stored. Default is 'Standard_Basis.h5'.

        Returns:
        --------
        StBas : np.ndarray
            The retrieved standard basis as a NumPy array.

        This function reads the standard basis from an HDF5 file and returns it as a NumPy array.

        """
        dim0 = 0
        dim1 = 0
        if (mpi.is_master_node()):
            ar = HDFArchive(h5name, 'r')
            if len(item) == 2:
                StBas = ar[item[0]][item[1]]
            else:
                StBas = ar[item[0]]
            dim0 = StBas.shape[0]
            dim1 = StBas.shape[1]
        mpi.barrier()
        dim0 = mpi.bcast(dim0)
        dim1 = mpi.bcast(dim1)
        if not mpi.is_master_node():
            StBas = np.zeros([dim0, dim1], np.complex)
        mpi.barrier()
        StBas = mpi.bcast(StBas)
        return StBas
