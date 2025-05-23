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
import triqs.utility.mpi as mpi
from h5 import *
from scipy import floor, sqrt
from scipy.special import factorial
import numpy as np

multipole_names = ['Monopole', 'Dipole', 'Quadrupole', 'Octupole', 'Hexadecapole'] + ['rank-%s' % (rank) for rank in
                                                                                      range(5, 16)]
real_harm = [['1'], ['y', 'z', 'x'], ['xy', 'yz', 'z2', 'xz', 'x2-y2'],
             ['y(x2-3y2)', 'xyz', 'yz2', 'z3', 'xz2', 'z(x2-y2)', 'x(3x2-y2)']]


class Multipolar:
    """This class manipulates multipole operators  ."""

    def __init__(self, J=None, Mult_Comp=None, action='dirprod', Kmax=None, basis='complex',
                 conv_mon_dipol='sph_tensor', print_mat=False, J1=None):
        """
        Initialize multipole tensors for magnetic interactions.

        This constructor initializes a Multipolar object, which can be used to represent magnetic interactions in various
        bases and multipole components. It allows for flexibility in specifying the properties of the multipole tensors.

        Parameters:
        -----------
        J : float, optional
            Total angular momentum quantum number. If provided, multipole tensors are generated based on this J value.
            Either `J` or `Mult_Comp` should be defined. Defaults to None.
        Mult_Comp : list of Multipolar objects, optional
            A list of two Multipolar objects. If provided, they are used to define a direct product of the multipole
            components. Either `J` or `Mult_Comp` should be defined. Defaults to None.
        action : str, optional
             The type of action to perform when initializing the multipole tensors. Can be 'dirprod' (direct product
             of two Multipolar objects) or 'merge'. Defaults to 'dirprod'.
        Kmax : int, optional
            Maximum K value for the multipole tensors. If not provided, it is calculated based on the value of `J`.
            Defaults to None.
        basis : str, optional
            The basis in which the multipole tensors are represented. Can be 'real' or 'complex'. Defaults to 'complex'.
        conv_mon_dipol : str, optional
            The convention used for monopoles and dipoles. Can be 'sph_tensor' or 'operator'.
            Defaults to 'sph_tensor'.
        print_mat : bool, optional
             If True, print information about the conversion of multipole tensors. Defaults to False.

        Notes:
        ------
        - The `action` parameter determines how the multipole tensors are initialized. For example, specifying 'dirprod'
          will initialize the object as a direct product of two `Mult_Comp` objects.
        - The multipole tensors are generated based on the specified `J` and `Kmax` values. If `J` is not provided, it
          can be calculated based on `Kmax`.
        """
        assert (J != None or Mult_Comp != None), 'Multipolar: either J or Mult_Comp should be defined!'

        def print_matrix(mat):
            ndim = mat.shape[0]
            ndim1 = mat.shape[1]
            for ind in range(ndim):
                str = ''
                for ind1 in range(ndim1):
                    str += "   %6.7f %6.7f" % (mat[ind, ind1].real, mat[ind, ind1].imag)
                print(str)

        def print_rmatrix(mat):
            ndim = mat.shape[0]
            ndim = mat.shape[1]
            for ind in range(ndim):
                str = ''
                for ind1 in range(ndim1):
                    str += " %4.1f" % (mat[ind, ind1])
                print(str)

        def print_tensors(J, Kmax, TKQ, conv_mon_dipol):
            print("\n\nTENSOR OPERATORS FOR J = %s" % (J))
            for K in range(Kmax + 1):
                print("\n\n%s operators:" % (multipole_names[K]))
                if K == 0 or K == 1:
                    if conv_mon_dipol == 'operator': print(
                        "Operators s(Q) are printed for %s" % (multipole_names[K]))
                for Q in range(-K, K + 1):
                    print("\nT(%s,%s): " % (K, Q))
                    mat = TKQ['%s_%s' % (K, Q)]
                    print_matrix(mat)

        def convert_to_real(Kmax, TKQ):
            for K in range(Kmax + 1):
                for Q in range(1, K + 1):
                    mmat = TKQ['%s_%s' % (K, Q)]
                    mmat1 = TKQ['%s_%s' % (K, -Q)]
                    mat = (mmat * np.power(-1, Q) + mmat1) * np.sqrt(0.5)
                    mat1 = (mmat - mmat1 * np.power(-1, Q)) * np.sqrt(-0.5 + 0j)
                    TKQ['%s_%s' % (K, Q)] = np.copy(mat)
                    TKQ['%s_%s' % (K, -Q)] = np.copy(mat1)

        def print_real_tensors(J, TKQ, Kmax):
            print("\n\nTENSOR OPERATORS in real basis FOR J = %s" % (J))
            for K in range(Kmax + 1):
                print("\n\n%s operators:" % (multipole_names[K]))
                if K == 0 or K == 1:
                    if conv_mon_dipol == 'moment': print(
                        "Moments J(Q) are printed for %s" % (multipole_names[K]))
                    if conv_mon_dipol == 'operator': print(
                        "Operators s(Q) are printed for %s" % (multipole_names[K]))
                for Q in range(-K, K + 1):
                    print("\nT(%s,%s): " % (K, Q))
                    if K < 4: print("Symmetry %s" % (real_harm[K][Q + K]))
                    mat = TKQ['%s_%s' % (K, Q)]
                    print_matrix(mat)

        if Mult_Comp == None:
            self.product = False
            self.conv_mon_dipol = conv_mon_dipol
            # generate from J
            self.J = J
            if not Kmax:
                self.Kmax = int(round(2 * J))
            else:
                self.Kmax = Kmax
            self.basis = basis

            if print_mat: print('conv_mon_dipol = %s' % conv_mon_dipol)
            TKQ = self.__gener_TKQ(J, conv_mon_dipol=conv_mon_dipol)
            if print_mat and mpi.is_master_node():
                print_tensors(J, self.Kmax, TKQ, conv_mon_dipol)

            if J1 != None:
                self.J1 = J1
                self.Kmax1 = int(round(2 * J1))
                if print_mat: print('conv_mon_dipol = %s' % conv_mon_dipol)
                TKQ1 = self.__gener_TKQ(J1, conv_mon_dipol=conv_mon_dipol)
                if print_mat and mpi.is_master_node():
                    print_tensors(J1, self.Kmax, TKQ, conv_mon_dipol)
            else:
                self.J1 = J
                self.Kmax1 = int(round(2 * J))
                TKQ1 = self.__gener_TKQ(J, conv_mon_dipol=conv_mon_dipol)

            if basis == 'real':
                convert_to_real(self.Kmax, TKQ)
                convert_to_real(self.Kmax1, TKQ1)
                if print_mat and mpi.is_master_node():
                    print_real_tensors(self.J, TKQ, self.Kmax)
                if J1 != None:
                    if print_mat and mpi.is_master_node():
                        print_real_tensors(self.J1, TKQ1, self.Kmax1)

            self.prn_bas = {}
            MKmax = max(self.Kmax, self.Kmax1)
            for K in range(MKmax + 1):
                for Q in range(-K, K + 1):
                    if K < 4 and basis == 'real':
                        self.prn_bas['%s_%s' % (K, Q)] = '{0: ^10}'.format(real_harm[K][K + Q])
                    else:
                        self.prn_bas['%s_%s' % (K, Q)] = '{0: ^10}'.format(Q)

            self.TKQ = TKQ
            self.TKQ1 = TKQ1

        elif action == 'dirprod':
            # Generate from the direct product
            self.product = True
            self.basis = Mult_Comp[0].basis
            Kmax0 = Mult_Comp[0].Kmax
            Kmax1 = Mult_Comp[1].Kmax
            J0 = Mult_Comp[0].J
            J1 = Mult_Comp[1].J
            TKQ = {}
            self.J = J0 * J1
            ndim = int(round(2 * self.J + 1))
            self.mult_names_new = {}
            self.Kmax = (Kmax0 + 1) * (Kmax1 + 1) - 1
            self.Qmax = {}
            if print_mat and mpi.is_master_node():
                print("\n\nTENSOR OPERATORS  FOR (Jxj)=(%s x %s)" % (J0, J1))
            for K0 in range(Kmax0 + 1):
                for K1 in range(Kmax1 + 1):
                    K = K0 * (Kmax1 + 1) + K1
                    self.mult_names_new[K] = self.__set_comp_label(K0, K1)
                    self.Qmax[K] = 2 * K0 * K1 + K0 + K1
                    for Q0 in range(-K0, K0 + 1):
                        for Q1 in range(-K1, K1 + 1):
                            Q = Q0 * (2 * K1 + 1) + Q1
                            mat0 = Mult_Comp[0].TKQ['%s_%s' % (K0, Q0)]
                            mat1 = Mult_Comp[1].TKQ['%s_%s' % (K1, Q1)]
                            TKQ['%s_%s' % (K, Q)] = np.kron(mat0, mat1)
                            if print_mat and mpi.is_master_node():
                                print("\nT(K,Q,k,q)=(%s,%s,%s,%s): " % (K0, Q0, K1, Q1))
                                print_matrix(TKQ['%s_%s' % (K, Q)])
            self.prn_bas = {}
            ind = 0
            for K0 in range(Kmax0 + 1):
                for Q0 in range(-K0, K0 + 1):
                    if K0 < 4 and Mult_Comp[0].basis == 'real':
                        prn = '[{0: ^3}'.format(real_harm[K0][K0 + Q0])
                    else:
                        prn = '[{0: ^3}'.format(Q0)
                    for K1 in range(Kmax1 + 1):
                        for Q1 in range(-K1, K1 + 1):
                            K = K0 * (Kmax1 + 1) + K1
                            Q = Q0 * (2 * K1 + 1) + Q1
                            if K1 < 4 and Mult_Comp[1].basis == 'real':
                                self.prn_bas['%s_%s' % (K, Q)] = prn + ';{0: ^3}] '.format(real_harm[K1][K1 + Q1])
                            else:
                                self.prn_bas['%s_%s' % (K, Q)] = prn + ';{0: ^3}] '.format(Q1)
            self.TKQ = TKQ
            self.ndim = ndim

        elif action == 'merge':

            self.product = True
            Kmax0 = Mult_Comp[0].Kmax
            Kmax1 = Mult_Comp[1].Kmax
            J0 = Mult_Comp[0].J
            J1 = Mult_Comp[1].J
            #assert Mult_Comp[0].stevens == Mult_Comp[
            #    1].stevens, "parameters stevens of shells J=%s and J=%s are not the same!" % (J0, J1)
            #self.stevens = Mult_Comp[0].stevens
            assert Mult_Comp[0].conv_mon_dipol == Mult_Comp[1].conv_mon_dipol, \
                "tensors for shells J=%s and J=%s generated with different conventions conv_mon_dipol!" % (J0, J1)
            self.conv_mon_dipol = Mult_Comp[0].conv_mon_dipol
            assert Mult_Comp[0].basis == Mult_Comp[1].basis, \
                "bases for shells J=%s and J=%s are different!" % (J0, J1)

            assert Mult_Comp[0].basis == 'real', \
                "both bases for shells J=%s and J=%s must be real!" % (J0, J1)
            self.basis = Mult_Comp[0].basis
            TKQ = {}
            ndim = int(round(2 * J0 + 2 * J1 + 2))
            ndim0 = int(round(2 * J0 + 1))
            ndim1 = int(round(2 * J1 + 1))
            self.mult_names_new = {}

            self.Kmax = Kmax0 + Kmax1 + 1
            self.Qmax = {}
            if print_mat and mpi.is_master_node():
                print("\n\nTENSOR OPERATORS in %s basis FOR (J\/j)=(%s \/ %s)" % (self.basis, J0, J1))
            self.prn_bas = {}
            for K0 in range(Kmax0 + 1):
                self.Qmax[K0] = K0
                self.mult_names_new[K0] = '[' + multipole_names[K0] + '(1)]'
                for Q in range(-K0, K0 + 1):
                    TKQ['%s_%s' % (K0, Q)] = np.zeros([ndim, ndim], complex)
                    TKQ['%s_%s' % (K0, Q)][0:ndim0, 0:ndim0] = Mult_Comp[0].TKQ['%s_%s' % (K0, Q)]
                    label = Mult_Comp[0].prn_bas['%s_%s' % (K0, Q)].rstrip().lstrip()
                    self.prn_bas['%s_%s' % (K0, Q)] = '[' + '{0: ^5}'.format(label) + '(1)]'
                    if print_mat and mpi.is_master_node():
                        print("\nT(K,Q)=(%s,%s), %s : " % (K0, Q, self.prn_bas['%s_%s' % (K0, Q)]))
                        print_matrix(TKQ['%s_%s' % (K0, Q)])
            for K1 in range(Kmax1 + 1):
                Ktot = Kmax0 + K1 + 1
                self.Qmax[Ktot] = K1
                self.mult_names_new[Ktot] = '[' + multipole_names[K1] + '(2)]'
                for Q in range(-K1, K1 + 1):
                    TKQ['%s_%s' % (Ktot, Q)] = np.zeros([ndim, ndim], complex)
                    TKQ['%s_%s' % (Ktot, Q)][ndim0:ndim, ndim0:ndim] = Mult_Comp[1].TKQ['%s_%s' % (K1, Q)]
                    label = Mult_Comp[1].prn_bas['%s_%s' % (K1, Q)].rstrip().lstrip()
                    self.prn_bas['%s_%s' % (Ktot, Q)] = '[' + '{0: ^5}'.format(label) + '(2)]'
                    if print_mat and mpi.is_master_node():
                        print("\nT(K,Q)=(%s,%s), %s : " % (Ktot, Q, self.prn_bas['%s_%s' % (Ktot, Q)]))
                        print_matrix(TKQ['%s_%s' % (Ktot, Q)])

            # add mixed multipoles
            TKQ_mixed_left = self.__gener_TKQ(J0, J1, conv_mon_dipol=self.conv_mon_dipol)
            TKQ_mixed_right = self.__gener_TKQ(J1, J0, conv_mon_dipol=self.conv_mon_dipol)
            if J0 > J1:
                prn_bas = Mult_Comp[0].prn_bas
            else:
                prn_bas = Mult_Comp[1].prn_bas

            Kmin_mix = abs(int(round(J0 - J1)))
            Kmax_mix = int(round(J0 + J1))
            sq2 = 1.0 / np.sqrt(2.0)
            isq2 = sq2 * 1j
            for K in range(Kmin_mix, Kmax_mix + 1):
                self.Kmax += 1
                self.Qmax[self.Kmax] = K
                # first kind: (T_KQ+T_KQ^dag)/sqrt(2)
                self.mult_names_new[self.Kmax] = 'mixed-positive ' + multipole_names[K]
                for Q in range(-K, K + 1):
                    TKQ['%s_%s' % (self.Kmax, Q)] = np.zeros([ndim, ndim], complex)
                    sign = 1 - 2 * (int(round(J0 - J1 + Q)) % 2)
                    TKQ['%s_%s' % (self.Kmax, Q)][ndim0:ndim, 0:ndim0] = TKQ_mixed_right['%s_%s' % (K, Q)] * sq2
                    TKQ['%s_%s' % (self.Kmax, Q)][0:ndim0, ndim0:ndim] = sign * TKQ_mixed_left['%s_%s' % (K, -Q)] * sq2
                    label = prn_bas['%s_%s' % (K, Q)].rstrip().lstrip() + '(+)'
                    self.prn_bas['%s_%s' % (self.Kmax, Q)] = '{0: ^10}'.format(label)
                # second kind: (T_KQ-T_KQ^dag)/*i/sqrt(2)
                self.Kmax += 1
                self.Qmax[self.Kmax] = K
                self.mult_names_new[self.Kmax] = 'mixed-negative ' + multipole_names[K]
                for Q in range(-K, K + 1):
                    TKQ['%s_%s' % (self.Kmax, Q)] = np.zeros([ndim, ndim], complex)
                    sign = 1 - 2 * (int(round(J0 - J1 + Q)) % 2)
                    TKQ['%s_%s' % (self.Kmax, Q)][ndim0:ndim, 0:ndim0] = TKQ_mixed_right['%s_%s' % (K, Q)] * isq2
                    TKQ['%s_%s' % (self.Kmax, Q)][0:ndim0, ndim0:ndim] = -sign * TKQ_mixed_left[
                        '%s_%s' % (K, -Q)] * isq2
                    label = prn_bas['%s_%s' % (K, Q)].rstrip().lstrip() + '(-)'
                    self.prn_bas['%s_%s' % (self.Kmax, Q)] = '{0: ^10}'.format(label)
                # combined moments to ensure a definite parity under time-reversal
                for Q in range(1, K + 1):
                    # Q> 0: first kind;  Q<0: second kind
                    mat = (TKQ['%s_%s' % (self.Kmax - 1, Q)] * np.power(-1, Q) + TKQ[
                        '%s_%s' % (self.Kmax - 1, -Q)]) * np.sqrt(0.5)
                    mat1 = (TKQ['%s_%s' % (self.Kmax, Q)] * np.power(-1, Q) - TKQ['%s_%s' % (self.Kmax, -Q)]) * np.sqrt(
                        0.5)
                    # vice versa
                    mat2 = (TKQ['%s_%s' % (self.Kmax, Q)] * np.power(-1, Q) + TKQ['%s_%s' % (self.Kmax, -Q)]) * np.sqrt(
                        0.5)
                    mat3 = (TKQ['%s_%s' % (self.Kmax - 1, Q)] * np.power(-1, Q) - TKQ[
                        '%s_%s' % (self.Kmax - 1, -Q)]) * np.sqrt(0.5)
                    TKQ['%s_%s' % (self.Kmax - 1, Q)] = np.copy(mat)
                    TKQ['%s_%s' % (self.Kmax - 1, -Q)] = np.copy(mat1)
                    TKQ['%s_%s' % (self.Kmax, Q)] = np.copy(mat2)
                    TKQ['%s_%s' % (self.Kmax, -Q)] = np.copy(mat3)
                # print
                if print_mat and mpi.is_master_node():
                    print("\n\n================== %ss:\n" % (self.mult_names_new[self.Kmax - 1]))
                    for Q in range(-K, K + 1):
                        print("\nT(K,Q)=(%s,%s), %s : " % (K0, Q, self.prn_bas['%s_%s' % (self.Kmax - 1, Q)]))
                        print_matrix(TKQ['%s_%s' % (self.Kmax - 1, Q)])
                    print("\n\n================== %ss:\n" % (self.mult_names_new[self.Kmax]))
                    for Q in range(-K, K + 1):
                        print("\nT(K,Q)=(%s,%s), %s : " % (K0, Q, self.prn_bas['%s_%s' % (self.Kmax, Q)]))
                        print_matrix(TKQ['%s_%s' % (self.Kmax, Q)])
            self.TKQ = TKQ
            self.TKQ1 = self.TKQ
            self.Kmax1 = self.Kmax
            self.J1 = J1

        else:
            print('Mult_Comp is given, but parameter action=%s is not recognized' % (action))
            print('action can be either \'dirprod\' or  \'merge\'')
            stop

        self.J1 = J1
        if self.J1 is None or action=='dirprod' :
            self.TKQ1 = self.TKQ
            self.Kmax1 = self.Kmax
            self.J1 = self.J

        self.TKQ_norm = {}
        self.TKQ1_norm = {}

        for key, mat in self.TKQ.items():
            if self.basis == 'real':
                self.TKQ_norm[key] = np.trace(np.dot(mat, mat))
            else:
                self.TKQ_norm[key] = np.trace(np.dot(mat, mat.conjugate().transpose()))
        for key, mat in self.TKQ1.items():
            if self.basis == 'real':
                self.TKQ1_norm[key] = np.trace(np.dot(mat, mat))
            else:
                self.TKQ1_norm[key] = np.trace(np.dot(mat, mat.conjugate().transpose()))

        if print_mat and mpi.is_master_node():
            print('\nNormalization of tensors:')
            for K in range(0, self.Kmax + 1):
                if self.product:
                    Qmax = self.Qmax[K]
                else:
                    Qmax = K
                for Q in range(-Qmax, Qmax + 1):
                    print('%s_%s, Norm = %s' % (K, Q, self.TKQ_norm['%s_%s' % (K, Q)]))
            if J1 != None:
                print('\nNormalization of tensors J1:')
                for K2 in range(0, self.Kmax1 + 1):
                    if self.product:
                        Qmax = self.Qmax[K2]
                    else:
                        Qmax = K2
                    for Q2 in range(-Qmax, Qmax + 1):
                        print('%s_%s, Norm = %s' % (K2, Q2, self.TKQ1_norm['%s_%s' % (K2, Q2)]))

    def Mult_interact(self, V4_ind, MagInt, fname='V_Mult.dat', tol_prnt=5e-7):
        """
        Transform V4 from atomic states to a multipolar basis and save the results to a file.

        This function transforms the interactions represented in atomic states (V4_ind) into a multipolar basis
        and stores the results in a file. It calculates the transformed interactions for various multipole components
        and coordinates, providing a detailed view of the transformed interactions.

        Parameters:
        -----------
        V4_ind : dict
             A dictionary containing atomic state interactions, typically generated from atomic calculations.
        MagInt : MagIntType
             An instance of the MagIntType class containing information about multipolar interactions.
        fname : str
             The name of the file where the transformed multipolar interactions will be saved. Default is 'V_Mult.dat'.
        tol_prnt : float
             A tolerance threshold for printing transformed interactions. Interactions with absolute values below this
             threshold will not be printed. Default is 1e-6.

        Returns:
        --------
        Vmult (dict): A dictionary containing the transformed multipolar interactions.
        """

        if (mpi.is_master_node()):
            f = open(fname, 'w')
            self.Vmult = {}
            for int in MagInt.Interactions:
                pair = MagInt.int_pairs[int]
                for key0 in pair.R_vec:
                    (type0, mul0), (type1, mul1) = key0
                    dim0 = MagInt.GSM[type0]
                    dim1 = MagInt.GSM[type1]
                    f.write("\n -----INTERACTION: %s_%s - %s_%s\n" % (type0, mul0, type1, mul1))
                    # for i in range(pair.nshells):
                    nsh_max = min(len(pair.coor_num[key0]), pair.nshells)
                    for i in range(nsh_max):
                        f.write("\nShell=%s\n" % (i))
                        f.write("\nCoor_num = %s\n" % (pair.coor_num[key0][i]))
                        # cycle over multipoles
                        for K in range(self.Kmax + 1):
                            if self.product:
                                Qmax = self.Qmax[K]
                            else:
                                Qmax = K
                            N = 2 * Qmax + 1
                            for K1 in range(self.Kmax1 + 1):
                                if self.product:
                                    Qmax1 = self.Qmax[K1]
                                else:
                                    Qmax1 = K1
                                N1 = 2 * Qmax1 + 1
                                str = '          '
                                for qq in range(-Qmax1, Qmax1 + 1): str += self.prn_bas['%s_%s' % (K1, qq)]
                                if self.product:
                                    f.write(
                                        "\n(%s)-(%s) interactions\n" % (
                                            self.mult_names_new[K], self.mult_names_new[K1]))
                                else:
                                    f.write("\n%s-%s interactions\n" % (multipole_names[K], multipole_names[K1]))
                                for i1 in range(pair.coor_num[key0][i]):
                                    R_key = "%s_%s" % (i, i1)
                                    key = "%s_%s_%s_%s_%s" % (type0, mul0, type1, mul1, R_key)
                                    if (key in V4_ind):
                                        V4 = V4_ind[key]
                                        m_key = key + "_%s_%s" % (K, K1)
                                        # f.write("m_key=%s\n"%(m_key))
                                        self.Vmult[m_key] = np.zeros([N, N1], np.complex)
                                        for Q in range(-Qmax, Qmax + 1):
                                            mat = self.TKQ['%s_%s' % (K, Q)] / self.TKQ_norm['%s_%s' % (K, Q)]
                                            for Q1 in range(-Qmax1, Qmax1 + 1):
                                                mat1 = self.TKQ1['%s_%s' % (K1, Q1)] / self.TKQ1_norm[
                                                    '%s_%s' % (K1, Q1)]
                                                for j in range(dim0):
                                                    for k in range(dim0):
                                                        for l in range(dim1):
                                                            for m in range(dim1):
                                                                self.Vmult[m_key][Qmax + Q, Qmax1 + Q1] += mat[j, k] * \
                                                                                                           V4[
                                                                                                               j, k, l, m] * \
                                                                                                           mat1[l, m]
                                        av_val = np.mean(abs(self.Vmult[m_key]))
                                        if av_val > tol_prnt:
                                            f.write("\n R = %s\n" % (pair.R_vec[key0][R_key]))
                                            f.write(" T = %s\n" % (np.array(pair.T_vec[key0][R_key])))
                                            # f.write(
                                            #    " D = %s\n" % (np.round(np.linalg.norm(pair.R_vec[key0][R_key])), 5))
                                            f.write("\n%s\n" % (str))
                                            prn_imag = False
                                            if np.mean(abs(self.Vmult[m_key].imag)) > tol_prnt: prn_imag = True
                                            for Q in range(-Qmax, Qmax + 1):
                                                st1 = "%s" % (self.prn_bas['%s_%s' % (K, Q)])
                                                for Q1 in range(-Qmax1, Qmax1 + 1):
                                                    if prn_imag:
                                                        st1 += '{0: ^11.4f}'.format(
                                                            1000.0 * self.Vmult[m_key][Qmax + Q, Qmax1 + Q1].real)
                                                        st1 += '{0: ^+8.4f}I'.format(
                                                            1000.0 * self.Vmult[m_key][Qmax + Q, Qmax1 + Q1].imag)
                                                    else:
                                                        st1 += '{0: ^10.4f}'.format(
                                                            1000.0 * self.Vmult[m_key][Qmax + Q, Qmax1 + Q1].real)
                                                f.write("%s\n" % (st1))
            f.close()

        return self.Vmult

    def store_results_in_h5(self, MagInt, fname):
        """
        Store essential data in an HDF5 file for future reference.

        This function is used to store important data related to multipolar interactions and other related information
        in an HDF5 file. The stored data can be used for future analysis, reference, or sharing with other users.

        Parameters:
        -----------
        MagInt : MagIntType
             An instance of the MagIntType class containing information about multipolar interactions.
        fname : str
             The name of the HDF5 file (excluding the '.h5' extension) where the data will be stored.
        """
        if (mpi.is_master_node()):
            ar = HDFArchive(fname + '.h5', 'a')
            if not ('Multipol_Int' in ar):
                ar.create_group('Multipol_Int')
            # data on types, interacting pairs, ground-state multiplicity
            ar['Multipol_Int']['interactions'] = MagInt.Interactions
            ar['Multipol_Int']['GSM'] = MagInt.GSM
            for int in MagInt.Interactions:
                int_key = '%s-%s' % (int[0], int[1])
                pairs = MagInt.int_pairs[int]
                ar['Multipol_Int']['shells_%s' % (int_key)] = pairs.nshells
                ar['Multipol_Int']['coor_nums_%s' % (int_key)] = pairs.coor_num
                ar['Multipol_Int']['R_vecs_%s' % (int_key)] = pairs.R_vec
                ar['Multipol_Int']['T_vecs_%s' % (int_key)] = pairs.T_vec
            # essential on tensors
            ar['Multipol_Int']['Kmax'] = self.Kmax
            ar['Multipol_Int']['ifproduct'] = self.product
            if self.product:
                ar['Multipol_Int']['Qmax'] = self.Qmax
                ar['Multipol_Int']['multipole_names'] = self.mult_names_new
            else:
                ar['Multipol_Int']['multipole_names'] = multipole_names
            ar['Multipol_Int']['prn_bas'] = self.prn_bas
            ar['Multipol_Int']['TKQ'] = self.TKQ
            # generalised
            if self.J != self.J1:
                ar['Multipol_Int']['TKQ1'] = self.TKQ1
                ar['Multipol_Int']['Kmax1'] = self.Kmax1
            # Vmult
            ar['Multipol_Int']['Vmult'] = self.Vmult
            del ar
            print("Multipolar interactions and other data stored in %s.h5" % (fname))

    def __gener_TKQ(self, J0, J1=None, conv_mon_dipol='sph_tensor'):
        """
        Generates operators matrices for given J0,J1.
        If only J0 is given, J1=J0 is assumed
        """
        TKQ = {}
        if J1 == None: J1 = J0
        ndim1 = int(round(2 * J0 + 1))
        ndim2 = int(round(2 * J1 + 1))
        Kmin = int(abs(round(J1 - J0)))
        Kmax = int(round(J1 + J0))
        for K in range(Kmin, Kmax + 1):
            for Q in range(-K, K + 1):
                TKQ['%s_%s' % (K, Q)] = np.zeros([ndim1, ndim2], complex)
                mat = TKQ['%s_%s' % (K, Q)]
                for ind in range(ndim1):
                    M = ind - J0
                    for ind1 in range(ndim2):
                        Mpr = ind1 - J1
                        mat[ind, ind1] += Wigner3j(J1, J0, K, Mpr, -M, Q) * np.power(-1, int(round(
                            J0 - M))) * np.sqrt(2.0 * K + 1.0)

        if (conv_mon_dipol == 'operator'):
            if (J0 != J1):
                print('Conv_mon_dipol =%s, but it must be =sph_tensor  for mixed multipoles (J0!=J1)!' \
                      % (conv_mon_dipol))
                assert 0
            # Convert TKQ to I for monopole
            TKQ['0_0'] *= np.sqrt(2 * J0 + 1)

            # Convert TKQ to s_Q operator for dipole moments
            if self.Kmax > 0:
                coff = np.sqrt(3.0 / ((2 * J0 + 1.0) * (J0 + 1.0) * J0))
                for Q in range(-1, 2):
                    TKQ['1_%s' % (Q)] /= coff

        elif (conv_mon_dipol != 'sph_tensor'):
            print('Multipolar: option %s for conv_mon_dipol not recognized!' % (conv_mon_dipol))
            assert 0
        return TKQ

    def __set_comp_label(self, K0, K1):
        """
        Combines two multipole names into a single one for composite multipole operators.
        This function is used to create a composite label for multipole operators, combining the labels of two
        multipole operators. The resulting label represents a composite operator formed by the interaction
        of two individual operators.

        Parameters:
        -----------
        K0 : int
             The label of the first multipole operator.
        K1 : int
             The label of the second multipole operator.

        Returns:
        --------
        comp_name (str):
             The composite label representing the interaction of the two multipole operators.

        """
        lK = [K0, K1]
        comp_name = ''
        for K in lK:
            if K == 2:
                str = multipole_names[K][0:4]
            elif K < 5:
                str = multipole_names[K][0:3]
            else:
                str = 'r-%s' % multipole_names[K][5:]
            comp_name += str
            if K == K0: comp_name += ';'
        return comp_name

    def __convert_to_stevens(self, TKQ, Kval=None):
        """
        Converts operators to Stevens operators using the normalization from Smith/Thornley (Proc. Phys. Soc. 89 p. 779)

        This function converts operators from a basis to Stevens operators, a commonly used basis set in quantum physics
        The normalization used is based on the work of Smith and Thornley, which ensures proper scaling.

        Parameters:
        -----------
        TKQ : dict
             A dictionary containing operators in the Tensorial Kramers-Quadratic (TKQ) notation.
        Kval : list, optional
             A list containing custom K values for each K index. If not provided, the function uses
        the standard K values from 0 to `self.Kmax`.

        Returns:
        --------
        None

        """
        for Kind in range(self.Kmax + 1):
            K = Kind if Kval == None else Kval[Kind]
            fact = np.sqrt(factorial(2 * self.J + K + 1) / factorial(2 * self.J - K)) / np.power(2,
                                                                                                 K) / np.sqrt(
                2 * K + 1.0)
            for Q in range(-K, K + 1):
                TKQ['%s_%s' % (Kind, Q)] *= fact

    def convert_to_mcphase_inp(self):
        """
        Convert Stevens operators in TKQ notation to the normalization used in the McPhase package.

        This function converts Stevens operators from TKQ (Tensorial Kramers-Quadratic) notation to the normalization
        used in the McPhase package. This conversion is necessary for consistency with McPhase's implementation.

        Parameters:
        -----------
        None

        Returns:
        --------
        None
        """
        a = [[1.0], [1.0, np.sqrt(0.5)], [0.5, np.sqrt(1.5) * 0.5, np.sqrt(3.0 / 8)],
             [0.5, np.sqrt(3.0 / 16) * 0.5, np.sqrt(15.0 / 8) * 0.5, np.sqrt(5.0 / 16)],
             [1.0 / 8, np.sqrt(5.0 / 16) * 0.5, np.sqrt(5.0 / 32) * 0.5, np.sqrt(35.0 / 16) * 0.5,
              np.sqrt(35.0 / 128)],
             [1.0 / 8, np.sqrt(15.0 / 128) * 0.5, np.sqrt(105.0 / 32) * 0.5, np.sqrt(35.0 / 256) * 0.5,
              np.sqrt(315.0 / 128) * 0.5, np.sqrt(63.0 / 256)],
             [1.0 / 16, np.sqrt(21.0 / 128) * 0.5, np.sqrt(105.0 / 1024) * 0.5, np.sqrt(105.0 / 256) * 0.5,
              np.sqrt(63.0 / 512) * 0.5, np.sqrt(693.0 / 256) * 0.5, np.sqrt(231.0 / 1024)]]
        b = [[1], [1, 0.5], [1, 0.25, 0.5], [1, 0.25, 0.25, 0.5], [1, 0.25, 0.25, 0.25, 0.5],
             [1, 0.25, 0.25, 0.25, 0.25, 0.5],
             [1, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5]]
        if self.Kmax > 6:
            Kconv = 6
            if mpi.is_master_node():
                print('Kmax = %s but conversion to mphase norms is done only up to Kmax=6' % (self.Kmax))
        else:
            Kconv = self.Kmax
        # for dipole mcphase I1,I2,I3 = T(1,1),T(1,-1),T(1,0) in Stevens not., hence skipped
        Karr = [0]
        Karr.extend(range(Kconv + 1))
        for K in Karr:
            self.TKQ['%s_0' % (K)] /= a[K][0]
            for Q in range(1, K + 1):
                cf = b[K][Q] / a[K][Q] * np.sqrt(2.0)
                print(K, Q, cf)
                self.TKQ['%s_%s' % (K, Q)] *= cf
                self.TKQ['%s_%s' % (K, -Q)] *= cf * np.power(-1.0, Q + 1)

    def upfold(self, trnmat):
        """
        Calculate upfolded interactions OT_kq O^herm_conj.

        This function calculates upfolded interactions between operators represented in the TKQ (Tensorial
        Kramers-Quadratic) notation. It is used to transform operators from the original basis to a new basis
        defined by the transformation matrix `trnmat`.

        Parameters:
        -----------
        trnmat (np.ndarray):
             The transformation matrix used to change the basis of the operators.

        Returns:
        --------
        None

        """
        assert trnmat.shape[1] == self.TKQ['0_0'].shape[0], 'upfold: dimension of trnmat not consistent with TKQ'
        ndim = trnmat.shape[0]
        for K in range(self.Kmax + 1):
            Qmax = K
            if self.product: Qmax = self.Qmax[K]
            for Q in range(-Qmax, Qmax + 1):
                mat_in = self.TKQ['%s_%s' % (K, Q)].copy()
                tmp = np.dot(trnmat, mat_in)
                self.TKQ['%s_%s' % (K, Q)] = np.dot(tmp, trnmat.transpose().conjugate())


def Wigner3j(j1, j2, j3, m1, m2, m3):
    """
    Compute the Wigner 3j symbol using the Racah formula.

    Parameters:
    -----------
    j1 : int or float
        The first total angular momentum quantum number.
    j2 : int or float
        The second total angular momentum quantum number.
    j3 : int or float
        The third total angular momentum quantum number.
    m1 : int or float
        The first magnetic quantum number.
    m2 : int or float
        The second magnetic quantum number.
    m3 : int or float
        The third magnetic quantum number.

    Returns:
    --------
    wigner3j : float
        The value of the Wigner 3j symbol for the specified quantum numbers.

    Usage:
    ------
    from wigner import Wigner3j
    wigner = Wigner3j(j1, j2, j3, m1, m2, m3)

    Symbol representation:
    ----------------------
     / j1  j2  j3 \
     |             |
     \ m1  m2  m3 /

    Notes:
    ------
    - The function uses the Racah formula for computation.
    - For more details, see the Wigner 3j-Symbol entry on Eric Weinstein's MathWorld:
      http://mathworld.wolfram.com/Wigner3j-Symbol.html
    """

    # Error checking
    if ((2 * j1 != floor(2 * j1)) | (2 * j2 != floor(2 * j2)) | (2 * j3 != floor(2 * j3)) | (
            2 * m1 != floor(2 * m1)) | (2 * m2 != floor(2 * m2)) | (2 * m3 != floor(2 * m3))):
        print('All arguments must be integers or half-integers.')
        return -1

    # Additional check if the sum of the second row equals zero
    if (m1 + m2 + m3 != 0):
        # print '3j-Symbol unphysical'
        return 0

    if (j1 - m1 != floor(j1 - m1)):
        # print '2*j1 and 2*m1 must have the same parity'
        return 0

    if (j2 - m2 != floor(j2 - m2)):
        # print '2*j2 and 2*m2 must have the same parity'
        return;
        0

    if (j3 - m3 != floor(j3 - m3)):
        # print '2*j3 and 2*m3 must have the same parity'
        return 0

    if (j3 > j1 + j2) | (j3 < abs(j1 - j2)):
        # print 'j3 is out of bounds.'
        return 0

    if abs(m1) > j1:
        # print 'm1 is out of bounds.'
        return 0

    if abs(m2) > j2:
        # print 'm2 is out of bounds.'
        return 0

    if abs(m3) > j3:
        # print 'm3 is out of bounds.'
        return 0

    t1 = j2 - m1 - j3
    t2 = j1 + m2 - j3
    t3 = j1 + j2 - j3
    t4 = j1 - m1
    t5 = j2 + m2

    tmin = max(0, max(t1, t2))
    tmax = min(t3, min(t4, t5))
    tvec = np.arange(tmin, tmax + 1, 1)

    wigner = 0

    for t in tvec:
        wigner += (-1) ** t / (factorial(t) * factorial(t - t1) * factorial(t - t2) * factorial(t3 - t) * factorial(
            t4 - t) * factorial(t5 - t))

    return wigner * (-1) ** (j1 - j2 - m3) * sqrt(
        factorial(j1 + j2 - j3) * factorial(j1 - j2 + j3) * factorial(-j1 + j2 + j3) / factorial(
            j1 + j2 + j3 + 1) * factorial(j1 + m1) * factorial(j1 - m1) * factorial(j2 + m2) * factorial(
            j2 - m2) * factorial(j3 + m3) * factorial(j3 - m3))
