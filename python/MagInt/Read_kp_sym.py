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
import triqs.utility.mpi as mpi
from numpy import *
import numpy as np
from triqs_dft_tools.converters.plovasp.vaspio import Kpoints
from triqs_dft_tools.converters.plovasp.vaspio import Poscar

class read_kp_sym:
    """
    Calculate the full Brillouin zone mesh from VASP and Wien2k.

    For Wien2k, the class reads k-points in the Irreducible Brillouin Zone (IBZ) from the `case.klist` file.
    The symmetry operations are taken from the `case.dmftsym` file to determine the full Brillouin zone mesh.

    For VASP, the full Brillouin zone is directly available in the `IBZKPT` file.
    """

    def __init__(self, general_par):
        """
        Initialize the class based on the specified DFT software.

        Parameters:
        -----------
            general_par (dict): Dictionary containing:

                - 'filename' (str): Name of the file.

                - 'dft_exec' (str): Name of the DFT software ('Wien2k' or 'Vasp').

                - 'folder' (str, optional): Relevant folder for VASP. Not used for Wien2k.

                - 'verbosity' (int): Verbosity level.

            filename (str, optional): File name, if specified. This argument is optional
                                      and defaults to None if not provided.
        """
        self.fname = general_par['filename']
        self.use_symmetries=general_par['use_symmetries']

        if general_par['dft_exec'] == 'Wien2k':
            self.__init_Wien2k_(general_par['verbosity'])
        elif general_par['dft_exec'] == 'Vasp':
            self.__init_Vasp_(general_par['folder'], general_par['verbosity'])
        # elif dft_exec == 'qe':
        #    __init_qe_(self, filename)

    def __init_Vasp_(self, pathfolder, verbosity):
        """
        Initialize data for VASP.

        This method is responsible for reading the IBZKPT file for VASP. It extracts the necessary
        information about the structure, k-points, and their associated properties to represent the Brillouin zone.

        Parameters:
        -----------
        pathfolder : str
            Path to the folder containing the VASP files.
        verbosity : int
            Level of verbosity during execution.

            - verbosity > 1: Reports the full Brillouin zone mesh details.

            - verbosity > 0: Reports the number of points in the full Brillouin zone.

        Attributes Set:
        ----------------
        kmesh_kpt : ndarray
            Array representing k-points.
        nktot : int
            Total number of k-points.
        K_shift : list of arrays
            Shift in k-points.
        syop_fbz : list
            Symmetry operations for full Brillouin zone.
        nk : int
            Number of k-points.
        k_fbz : list of arrays
            k-points in the full Brillouin zone.
        flat_k_ind : list of lists
            Flat index for k-points.
        weight : float
            Weight associated with each k-point.

        Notes:
        ------
        The method also uses 'Kpoints' class.
        """

        filename = 'IBZKPT'
        KP = Kpoints()
        KP.from_file(vasp_dir=pathfolder, ibz_filename=filename)
        self.KP = KP

        self.kmesh_kpt = 0.0 * KP.kpts
        self.nktot = KP.nktot

        for ik in range(self.nktot):
            # for ii in range(3):
            self.kmesh_kpt[ik] += KP.kpts[ik].T  # @ (self.struct.kpt_basis[:] @ KP.kpts[ik]).T).T

        self.K_shift = [[np.array(x) for x in np.zeros((KP.nktot, 3), dtype=float).tolist()]]
        self.syop_fbz = [[i for i in range(KP.nktot)]]
        self.nk = 0
        # self.k_fbz = KP.kpts
        # self.k_fbz = self.k_fbz.tolist()
        # self.k_fbz = [[np.array(x) for x in self.k_fbz]]

        self.k_fbz = [[np.array(x) for x in self.kmesh_kpt.tolist()]]

        if (verbosity > 1):
            mpi.report("\nFull BZ mesh:")
            mpi.report("\n       #   #-IBZ              k-point              N_sym              K-shift")
        self.flat_k_ind = []
        for i in range(len(self.k_fbz)):
            self.flat_k_ind.append([])
            # for jj in range(3):
            kp_i = self.k_fbz[i]
            K_i = self.K_shift[i]
            syop_i = np.zeros(self.nktot)
            for ik in range(len(kp_i)):
                self.flat_k_ind[-1].append(self.nk)
                self.nk += 1
                if (verbosity > 1):
                    x = kp_i[ik][0]
                    y = kp_i[ik][1]
                    z = kp_i[ik][2]
                    mpi.report("  %6i  %4i       %8.5f %8.5f %8.5f  %5i      %8.5f %8.5f %8.5f" % (
                        self.nk, ik, x, y, z, syop_i[ik], K_i[ik][0], K_i[ik][1], K_i[ik][2]))

        self.weight = 1.0 / self.nk
        if verbosity > 0:
            mpi.report("\nNumber of points in the full BZ = %s" % (self.nk))

    def __init_Wien2k_(self, verbosity=0):
        """
        Initializes the class with data from Wien2k output files. This method reads symmetry operations, transformation
        matrices, and k-vectors, and then expands the k-point mesh to the full Brillouin zone (BZ).

        Parameters
        ----------
        verbosity : int, optional
            Level of verbosity for reporting. Default is 0, which means no reporting.

        Returns
        -------
        None

        Notes
        -----
        - Reads symmetry operations from a ".dmftsym" file.
        - Reads the transformation matrix from internal to Cartesian coordinates from the "outputkgen" file.
        - Reads k-vectors in internal coordinates from the "outputkgen" file.
        - Expands the mesh to the full Brillouin zone using the `__k_star_` function.

        Raises
        ------
        IOError
            If any of the required files cannot be opened.

        """
        self.symop = []
        self.symprop = []
        if self.use_symmetries:
            f = open(self.fname + ".dmftsym")
            isym = 0
            mat = zeros([3, 3])
            for il in f:
                if il.find('Global') > -1:
                    break
                if il.find('euler') > -1:
                    self.symprop.append(int(il.split()[3]))
                    irow = 0
                    for il1 in f:
                        mat[irow, 0:3] = [float(il1.split()[i]) for i in range(0, 3)]
                        irow += 1
                        if irow == 3:
                            break
                    self.symop.append(mat.copy())
            f.close()
        else:
            self.symop.append(np.identity(3,complex))
            self.symprop.append(1)

        # Reads the transformation from inner to cartesian coordinates from 'outputkgen'
        f = open(self.fname + ".outputkgen")
        self.G_trans = zeros([3, 3])
        for il in f:
            if il.find('G1        G2        G3') > -1:
                irow = 0
                for il1 in f:
                    self.G_trans[irow, 0:3] = [float(il1.split()[i]) for i in range(0, 3)]
                    irow += 1
                    if irow == 3:
                        break
                break
        f.close()
        # Reads k-vectors in internal coordinates
        f = open(self.fname + ".outputkgen")
        ik = 0
        kvec = []
        kvec_car = []
        vec = zeros([3])
        self.nk = 0
        for il in f:
            if il.find('internal and cartesian k-vectors:') > -1:
                try:
                    for il1 in f:
                        vec = [float(il1.split()[i]) for i in range(3)]
                        kvec.append(vec)
                        ik += 1
                        kvec_car.append(dot(self.G_trans, vec))
                except:
                    self.nk = ik
                    break
        self.G_trans_inv = linalg.inv(self.G_trans)
        self.kvek_ibz = kvec_car
        f.close()
        if verbosity > 0:
            mpi.report("\nNumber of points in the IBZ = %s\n" % (self.nk))
            mpi.report("Internal and cartesian coordinates:")
            for i in range(len(kvec)):
                k1 = kvec[i]
                k2 = kvec_car[i]
                mpi.report(
                    "%6i %8.5f %8.5f %8.5f       %8.5f %8.5f %8.5f" % (i, k1[0], k1[1], k1[2], k2[0], k2[1], k2[2]))
        # Expand the mesh to full BZ
        self.k_fbz = []
        self.syop_fbz = []
        self.K_shift = []
        for kp in range(len(kvec)):
            t_k, t_syop, t_shift = self.__k_star_(kvec[kp], kvec_car[kp])
            self.k_fbz.append(t_k)
            self.syop_fbz.append(t_syop)
            self.K_shift.append(t_shift)
        self.nk = 0
        if (verbosity > 0):
            mpi.report("\nFull BZ mesh:")
            mpi.report("\n       #   #-IBZ              k-point              N_sym              K-shift")
        self.flat_k_ind = []
        for i in range(len(self.k_fbz)):
            self.flat_k_ind.append([])
            kp_i = self.k_fbz[i]
            K_i = self.K_shift[i]
            syop_i = self.syop_fbz[i]
            # mpi.report(" k-star: %s"%(i))
            for ik in range(len(kp_i)):
                self.flat_k_ind[-1].append(self.nk)
                self.nk += 1
                # DEBUG
                # kp_i[ik] = kp_i[ik]
                # END DEBUG
                kp_i[ik] = dot(self.G_trans, kp_i[ik])
                K_i[ik] = dot(self.G_trans, K_i[ik])
                if verbosity > 0:
                    x = kp_i[ik][0]
                    y = kp_i[ik][1]
                    z = kp_i[ik][2]
                    mpi.report("  %6i  %4i       %8.5f %8.5f %8.5f  %5i      %8.5f %8.5f %8.5f" % (
                        self.nk, ik, x, y, z, syop_i[ik], K_i[ik][0], K_i[ik][1], K_i[ik][2]))

        self.weight = 1.0 / self.nk
        if verbosity > 0:
            mpi.report("\nNumber of points in the full BZ = %s" % (self.nk))

    def __k_star_(self, kvec, kvec_car, tol=1e-4):
        """
        Generates a k-points star from a given k-point in the irreducible Brillouin zone using symmetry operations.

        Parameters
        ----------
        kvec : array-like
            The k-point in fractional coordinates within the IBZ.
        kvec_car : array-like
            The k-point in Cartesian coordinates.
        tol : float, optional
            Tolerance for numerical precision in determining equivalent k-points. Default is 1e-4.

        Returns
        -------
        t_k : list of array-like
            List of k-points in fractional coordinates that form the star.
        t_syop : list of int
            List of symmetry operation indices corresponding to each k-point in the star.
        t_shift : list of array-like
            List of shift vectors applied to each k-point to bring it back into the first Brillouin zone.

        Raises
        ------
        AssertionError
            If a duplicate k-point is found within the tolerance `tol`.
        """
        t_k = []
        t_syop = []
        t_shift = []
        for isym in range(len(self.symop)):
            mat = self.symop[isym]
            knew_car = dot(mat, kvec_car) * self.symprop[isym]
            knew = dot(self.G_trans_inv, knew_car)
            K = array([0.0, 0.0, 0.0])
            try:
                for i in range(3):
                    while (knew[i] > 1.0 - tol):
                        knew[i] -= 1.0
                        K[i] += 1.0
                    while (knew[i] < 0.0 - tol):
                        knew[i] += 1.0
                        K[i] -= 1.0
                for vec in t_k:
                    dd = 0.0
                    for i in range(3):
                        dd += abs(vec[i] - knew[i])
                    if (dd < tol):
                        assert 0
            except:
                continue
            t_k.append(knew)
            t_syop.append(isym)
            t_shift.append(K)
        return t_k, t_syop, t_shift
