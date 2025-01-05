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
import numpy as np
import triqs.utility.mpi as mpi
import math
import re


def print_arr(arr, log='arr', fname=None, decd=4, prn_zero_imag=False):
    """
    Print the real and/or imaginary parts of an array with customizable format.

    This function is designed to print either the real part, imaginary part,
    or both of an array, depending on the presence of non-zero values in these
    parts and the parameters provided.

    Parameters
    ----------
    arr : numpy.ndarray
        The array to be printed. Should be a two-dimensional array.
    log : str, optional
        Label for the output. Defaults to 'arr'.
    fname : str, optional
        If provided, the output will be written to a file with this name.
        Otherwise, it will be printed to the standard output. Defaults to None.
    decd : int, optional
        Number of decimal places for the printed values. Defaults to 4.
    prn_zero_imag : bool, optional
        If True, the imaginary part of the array will always be printed,
        even if it's zero. Defaults to False.
    """
    m, n = arr.shape
    imag = False
    dig = decd + 5
    if not fname:
        print(log)
    else:
        f = open(fname, 'w')
        if log != 'arr': f.write("# %s \n" % log)
    if not prn_zero_imag:
        cmplx = False
        if (np.sum(abs(arr.imag)) > 1e-8):
            if (np.sum(abs(arr.real)) < 1e-8):
                imag = True
                if mpi.is_master_node():
                    print('Imaginary part:')
            else:
                if mpi.is_master_node():
                    cmplx = True
                    print('Complex:')
        else:
            if mpi.is_master_node():
                print('Real part:')
            imag = False
    else:
        cmplx = True
    for i in range(m):
        str = ''
        for j in range(n):
            if cmplx:
                str += ' %*.*f %*.*f' % (dig, decd, arr[i, j].real, dig, decd, arr[i, j].imag)
            elif imag:
                str += ' %*.*f' % (dig, decd, arr[i, j].imag)
            else:
                str += ' %*.*f' % (dig, decd, arr[i, j].real)
        if not fname:
            if mpi.is_master_node():
                print(str)
        else:
            f.write("%s\n" % str)
    if fname: f.close()


def print_mat(mat):
    """
    Print the real and imaginary parts of a matrix in a formatted manner.

    This function is designed to display each element of the matrix with its
    real and imaginary components. The output is formatted to show up to 5 decimal
    places for both the real and imaginary parts.

    Parameters
    ----------
    mat : numpy.ndarray
        The matrix to be printed. Should be a square matrix with both real and imaginary parts.

    Notes
    -----
    The output format for each element is: "real_part imaginary_part", with each part displayed with up to 5 decimal
    places.
    """
    ndim = np.shape(mat)[0]
    for ind in range(ndim):
        str = ''
        for ind1 in range(ndim):
            str += "%6.5f %6.5f " % (mat[ind, ind1].real, mat[ind, ind1].imag)
        print(str)


def read_eal(eal_fname, nlms):
    """
    Read the Effective Atomic Levels (EAL) from a file and store it in a dictionary format.

    This function reads the EAL data from the specified file. The EAL data is expected
    to contain real and imaginary parts, which are combined into a complex number.
    The processed data is then stored in a dictionary with a key 'ud'.

    Parameters
    ----------
    eal_fname : str
        Path to the file containing the EAL data.

    nlms : int
        Dimension of the EAL square matrix. It indicates the number of rows and columns of the matrix.

    Returns
    -------
    dict
        Dictionary containing the EAL data as a complex matrix under the 'ud' key.
    """
    # Load Effective Atomic Levels
    eal_tmp = np.loadtxt(eal_fname)
    eal = {'ud': np.zeros((nlms, nlms), complex)}
    for i in range(nlms):
        for j in range(nlms):
            eal['ud'][i, j] = eal_tmp[i, 2 * j] + 1j * eal_tmp[i, 2 * j + 1]
    return eal


def read_mat(mat_fname):
    """
    Read a matrix from a file and convert it to a complex format.
 
    This function reads data from the specified file where the data is expected to have real and
    imaginary parts presented in separate columns. The function combines these parts into complex numbers
    and arranges them into a complex matrix of appropriate dimensions.
 
    Parameters
    ----------
    mat_fname : str
        Path to the file containing the matrix data.
 
    Returns
    -------
    np.ndarray
        A complex matrix constructed from the real and imaginary parts read from the file.
 
    """

    # Load Effective Atomic Levels
    mat_tmp = np.loadtxt(mat_fname)
    dim1 = int(mat_tmp.shape[0])
    dim2 = int(mat_tmp.shape[1] / 2)
    mat = np.zeros((dim1, dim2), complex)
    for i in range(dim1):
        for j in range(dim2):
            mat[i, j] = mat_tmp[i, 2 * j] + 1j * mat_tmp[i, 2 * j + 1]
    return mat


def gener_S_mat(nlm):
    """
    Generate the spin matrices Sx, Sy, Sz, Sp and Sm, where Sp/Sm are the S+/- ladder operators respectively.

    This function creates the three spin matrices for a given dimension.
    The matrices are constructed in a larger 2D space, where the original dimension is multiplied by 2
    to account for spin up and spin down states. Each spin matrix is represented in the Pauli matrix form.

    Parameters
    ----------
    nlm : int
        Dimension for which the spin matrices are to be constructed.

    Returns
    -------
    dict
        A dictionary containing the spin matrices 'x', 'y', 'z', 'p', 'm' , where each key maps to its respective
        spin matrix.
    """
    S = {}
    S['z'] = np.zeros((2 * nlm, 2 * nlm), complex)
    S['p'] = np.zeros((2 * nlm, 2 * nlm), complex)
    S['m'] = np.zeros((2 * nlm, 2 * nlm), complex)
    for i in range(nlm):
        S['p'][i, i + nlm] = 1.0
        S['m'][i + nlm, i] = 1.0
        S['z'][i, i] = 0.5
        S['z'][i + nlm, i + nlm] = -0.5
    S['x'] = 0.5 * (S['p'] + S['m'])
    S['y'] = 0.5j * (S['m'] - S['p'])
    return S


def gener_J_mat(j):
    Ndim = int(2 * j) + 1
    J = {}
    J['z'] = np.zeros((Ndim, Ndim), complex)
    J['p'] = np.zeros((Ndim, Ndim), complex)
    J['m'] = np.zeros((Ndim, Ndim), complex)
    for i in range(Ndim):
        # m=j-i
        m = i - j
        J['z'][i, i] = m
        if i > 0:
            J['p'][i, i - 1] = math.sqrt((j - m + 1.0) * (j + m))
            J['m'][i - 1, i] = math.sqrt((j - m + 1.0) * (j + m))
    return J


def gener_L_mat(l):
    """
    Generate the total angular momentum matrices Jx, Jy, Jz, J+, and J-.

    This function constructs the total angular momentum matrices (Jx, Jy, Jz, J+, and J-)
    for a given quantum number j. The matrices are based on the quantum mechanics definitions
    for angular momentum.

    Parameters
    ----------
    j : float
        The total angular momentum quantum number.

    Returns
    -------
    dict
        A dictionary containing the angular momentum matrices 'x', 'y', 'z', 'p' (which corresponds to J+), and 'm'
        (which corresponds to J-), where each key maps to its respective matrix.
    """
    ll = gener_J_mat(l)
    L = {}
    L['z'] = np.kron(np.identity(2), ll['z'])
    L['p'] = np.kron(np.identity(2), ll['p'])
    L['m'] = np.kron(np.identity(2), ll['m'])
    L['x'] = 0.5 * (L['p'] + L['m'])
    L['y'] = 0.5j * (L['m'] - L['p'])
    return L


def gener_M_mat(L, S, gamma):
    """
    Generate magnetic moment operator matrices.

    Constructs the magnetic moment operators (Mx, My, and Mz) from the given
    orbital (L) and spin (S) angular momentum matrices using the Landé g-factor.

    Parameters
    ----------
    L : dict
        Dictionary containing the orbital angular momentum matrices 'x', 'y', and 'z'.

    S : dict
        Dictionary containing the spin angular momentum matrices 'x', 'y', and 'z'.

    gamma : float
        Landé g-factor used to weight the contribution from orbital angular momentum.

    Returns
    -------
    dict
        A dictionary containing the magnetic moment matrices 'x', 'y', and 'z'.
    """
    M = {}
    M['x'] = gamma * L['x'] + 2.0 * S['x']
    M['y'] = gamma * L['y'] + 2.0 * S['y']
    M['z'] = gamma * L['z'] + 2.0 * S['z']
    return M


def my_get_vec(n, iv, vec):
    """
    Extracts a specific vector from the 'STATES' file.

    Reads the 'STATES' file and retrieves the vector number `iv`.
    The structure of each line in 'STATES' is expected to have the format
    of a complex number in parenthesis, e.g., `(real,imag)`.

    Parameters
    ----------
    n : int
        The total number of vectors in the 'STATES' file.

    iv : int
        The index of the vector to be extracted (0-based index).

    vec : list or ndarray
        An initialized list or array to store the extracted vector values.

    Returns
    -------
    list or ndarray
        The updated list or array containing the values of the extracted vector.

    """
    # Read and return vector number iv from file 'STATES'
    with open('STATES', 'r') as f:
        for i in range(n):
            for k in range(n):
                line = f.readline()
                line = re.split(r'[(,)]', line)
                real = float(line[1])
                imag = float(line[2])
                vec[k] = complex(real, imag)
            #     print(vec[k])
            # print('')
            if i == iv:
                break
    return vec


def cub_eal(dim):
    """
    Extracts the transformation matrix from the 'real_d_harms' file.

    Reads the 'real_d_harms' file and constructs the transformation matrix.

    Parameters
    ----------
    dim : int
        The dimension of the square transformation matrix to be constructed.

    Returns
    -------
    ndarray
        The transformation matrix of shape (dim, dim) populated with complex values.
    """
    tt = np.loadtxt('real_d_harms')
    transmat = np.zeros((dim, dim), complex)
    for i in range(dim):
        for j in range(dim):
            transmat[i, j] = tt[i, 2 * j] + 1j * tt[i, 2 * j + 1]
    return transmat


def generate_PS_mat(S, T, orbital, psBas=None, ps=True):
    if ps == True:

        T['z'] = np.dot(psBas.conjugate().transpose(), np.dot(T['z'], psBas))
        T['p'] = np.dot(psBas.conjugate().transpose(), np.dot(T['p'], psBas))
        T['m'] = np.dot(psBas.conjugate().transpose(), np.dot(T['m'], psBas))

        lz = np.kron(np.identity(orbital), T['z'])
        lp = np.kron(np.identity(orbital), T['p'])
        lm = np.kron(np.identity(orbital), T['m'])

        J = {}
        J['z'] = S['z'] + lz
        J['p'] = S['p'] + lp
        J['m'] = S['m'] + lm
        J['x'] = 0.5 * (J['p'] + J['m'])
        J['y'] = 0.5j * (J['m'] - J['p'])
        return J, lz
    else:
        J = {}
        J['z'] = S['z'] + T['z']
        J['p'] = S['p'] + T['p']
        J['m'] = S['m'] + T['m']
        J['x'] = 0.5 * (J['p'] + J['m'])
        J['y'] = 0.5j * (J['m'] - J['p'])
        return J


def Kanamori_hamiltonian(U_int, J_hund, orbital, dft_exec):
    """
    Constructs the Kanamori Hamiltonian in the Wien2k format or returns the double counting
    for the Kanamori convention.

    Parameters
    ----------
    U_int : float
        The screened Coulomb interaction value.

    J_hund : float
        The Hund's exchange coupling value.

    orbital : int
        Defines the orbital angular momentum quantum number (l).
        For instance, for d orbitals, this would be 2.

    dft_exec : str
        The DFT executable being used. Currently, only 'Wien2k' is supported.

    Returns
    -------
    U_full : ndarray
        The full Kanamori Hamiltonian matrix in the given DFT executable's format.

    U_dc : float
        The double counting term for the Coulomb interaction.

    J_dc : float
        The double counting term for Hund's exchange interaction.
    """
    UK = U_int + 8 * J_hund / 7.0
    JK = 0.772 * J_hund
    JP = 1.0 * JK

    U_dc = UK - 2 * (JK + JP) / 3.0
    J_dc = 7 * JK / 3.0 - 2 * JP / 3.0

    nlm = 2 * orbital + 1
    nlms = nlm * 2

    # set matrix elements of Kanamory U
    U_Kan = np.zeros((nlm, nlm, nlm, nlm))
    U_full = np.zeros((nlms, nlms, nlms, nlms))

    if dft_exec == 'Wien2k':
        for i in range(nlm):
            for j in range(nlm):
                if i == j:
                    U_Kan[i, i, i, i] = UK
                else:
                    U_Kan[i, j, i, j] = UK - JK - JP
                    U_Kan[i, j, j, i] = JK
                    if i == 1 or j == 1:
                        U_Kan[i, i, j, j] = -JP
                    else:
                        U_Kan[i, i, j, j] = JP

        for sp1 in range(2):
            for sp2 in range(2):
                U_full[sp1 * nlm:(sp1 + 1) * nlm, sp2 * nlm:(sp2 + 1) * nlm, sp1 * nlm:(sp1 + 1) * nlm,
                sp2 * nlm:(sp2 + 1) * nlm] = U_Kan

    return U_full, U_dc, J_dc


def wien2k_basis(options='t2g'):
    if options == 'full_d':
        basis = read_mat('read_harm/full_d')
    elif options == 't2g':
        basis = read_mat('read_harm/real_p_harms')
        basis = -basis
    elif options == 'p':
        basis = read_mat('read_harm/real_p_harms')
        # basis[0, :] = 1j * basis[0, :]
        basis[1, :] = 1j * basis[1, :]
        # basis[3, :] = 1j * basis[3, :]
        basis[4, :] = 1j * basis[4, :]
    elif options == 'd':
        basis = read_mat('read_harm/real_d_harms')
    else:
        basis = read_mat('read_harm/real_f_harms')
    return basis


def compare_lists_arrays(list1, list2, rtol=1e-5, atol=1e-8):
    """
    Compare two lists of NumPy arrays.

    Parameters
    ----------

    list1: list
        First list of NumPy arrays.
    list2: list
        Second list of NumPy arrays.
    rtol: float
        Relative tolerance.
    atol: float
        Absolute tolerance.

    Returns
    -------
    True if all corresponding arrays in the lists are equal within the given tolerance, False otherwise.
    """

    # Check if lists have the same length
    if len(list1) != len(list2):
        return False

    # Compare arrays in the lists
    for arr1, arr2 in zip(list1, list2):
        if not np.allclose(arr1, arr2, rtol=rtol, atol=atol):
            return False

    return True


def compare_dicts(dict1, dict2, rtol=1e-5, atol=1e-8):
    """
        Compare two dictionaries for equality, with support for numpy arrays.

        Parameters
        ----------
        dict1 : dict
            The first dictionary to compare.
        dict2 : dict
            The second dictionary to compare.
        rtol : float, optional
            The relative tolerance parameter for numpy array comparison (default is 1e-5).
        atol : float, optional
            The absolute tolerance parameter for numpy array comparison (default is 1e-8).

        Returns
        -------
        bool
            Returns True if both dictionaries are the same; otherwise, False.

        """
    if dict1.keys() != dict2.keys():
        return False
    for key in dict1:
        val1 = dict1[key]
        val2 = dict2[key]
        if isinstance(val1, list) and all(isinstance(i, np.ndarray) for i in val1):
            if not compare_lists_of_arrays(val1, val2, rtol, atol):
                return False
        elif isinstance(val1, np.ndarray) and isinstance(val2, np.ndarray):
            if not np.array_equal(val1, val2):
                return False
        elif val1 != val2:
            return False
    return True


def compare_lists(list1, list2, rtol=1e-5, atol=1e-8):
    return len(list1) == len(list2) and all(np.isclose(a, b, rtol=rtol, atol=atol) for a, b in zip(list1, list2))


def is_unitary(matrix, tolerance=1e-4):
    identity = np.eye(matrix.shape[0])
    product = np.dot(matrix, np.conjugate(matrix.T))
    print("\nDeviation from unitary")
    print_arr(identity - np.dot(matrix, np.conjugate(matrix.T)))
    return np.allclose(product, identity, atol=tolerance)


def print_logo():
    mpi.report(' __  __   __   __  __  _  _  ____')
    mpi.report('(  \/  ) (  ) / _)(  )( \( )(_  _)')
    mpi.report(' )    (  /__\( (/\ )(  )  (   )(  ')
    mpi.report('(_/\/\_)(_)(_)\__/(__)(_)\_) (__) \n')

def print_warning_vasp():
    mpi.report('\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    mpi.report('!!     ▗▖ ▗▖ ▗▄▖ ▗▄▄▖ ▗▖  ▗▖▗▄▄▄▖▗▖  ▗▖ ▗▄▄▖     !!')
    mpi.report('!!     ▐▌ ▐▌▐▌ ▐▌▐▌ ▐▌▐▛▚▖▐▌  █  ▐▛▚▖▐▌▐▌        !!')
    mpi.report('!!     ▐▌ ▐▌▐▛▀▜▌▐▛▀▚▖▐▌ ▝▜▌  █  ▐▌ ▝▜▌▐▌▝▜▌     !!')
    mpi.report('!!     ▐▙█▟▌▐▌ ▐▌▐▌ ▐▌▐▌  ▐▌▗▄█▄▖▐▌  ▐▌▝▚▄▞      !!')
    mpi.report('!!                                               !!')
    mpi.report('!!  The VASP interface is experimental and has   !!')
    mpi.report('!!  not been extensively tested.                 !!')
    mpi.report('!!  Please contact the authors to get assistance !!')
    mpi.report('!!                                               !!')
    mpi.report('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
