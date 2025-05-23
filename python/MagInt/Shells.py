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
from numpy import *
import triqs.utility.mpi as mpi
from triqs_dft_tools.converters.plovasp.vaspio import Poscar
from MagInt.utils import *
from MagInt.Parser import *
import math
import os
import sys
import numpy as np
import itertools


class shells:
    """
    Initialize the class for analyzing interacting shells in various DFT software outputs.
    """

    def __init__(self, general_par, filename=None, interacting_types=None, interact_sites=None,
                 interact_basis=None, max_n_shells=3, print_vecs='Translation', R_arr=None, R_arr_site=None,
                 R_arr_origin=None, run_only=False, calc_all_sites=False):
        """
        Initialize the class for analyzing interacting shells in various DFT software outputs.

        Parameters:
        -----------
        general_par : dict
            General parameters dictionary containing:

            - 'dft_exec' (str): Specifies the DFT software ('Wien2k', 'Vasp', etc.).

            - 'verbosity' (int): Specifies the verbosity level.

            - 'tol_shells' (int): tollerance N for vectors length difference 10^(-N)
                                      in sorting them between coordination shells

        filename : str, optional
            File name for the data. Defaults to None.
        interacting_types : list, optional
            List of two atom types that are interacting. Defaults to None.
        interact_sites : list, optional
            List of interacting sites. Defaults to None.
        interact_basis : list, optional
            Basis for the interaction. Defaults to None.
        max_n_shells : int, optional
            Maximum number of shells to consider. Defaults to 3.
        R_arr : list, optional
            Custom array of interaction vectors. Defaults to None.
        R_arr_site : list, optional
            Sites associated with R_arr. Defaults to None.
        R_arr_origin : list, optional
            Origin points for R_arr. Defaults to None.
        calc_all_sites : bool
            True: lattice vectors connecting all sites of given types are calculated
            False : only lattice vectors connecting the representative sites with all other sites are calculated

        Notes:
        ------
        This method initializes the class based on the DFT software specified in `general_par['dft_exec']`.
        Currently supported DFT software are 'Wien2k' and 'Vasp'. Support for 'QuantumEspresso' may be added in the future.
        """

        self.atom0 = interacting_types[0]
        self.atom1 = interacting_types[1]
        self.single_atom = True
        if interacting_types[0] != interacting_types[1]:
            self.single_atom = False
        self.interact_sites = interact_sites
        self.interact_basis = interact_basis
        self.R_vec = {}
        self.R_site = {}
        self.R_shell = []
        self.coor_num = []
        self.R_origin = {}
        self.R_basis = {}
        self.R_arr_set = False
        self.atom_list = None
        self.name_atoms = None
        self.calc_all_sites = calc_all_sites
        direct = False
        tol_shells = general_par['tol_shells']

        # Initialize Wien2k shells
        if general_par['dft_exec'] == 'Wien2k':
            self.__init_Wien2k_(filename)
            direct = True

        # Initialize VASP shells
        elif general_par['dft_exec'] == 'Vasp':
            self.__init_Vasp_(general_par)
            direct = True

        # Initialize Quantum Espresso shells
        # elif general_par['dft_exec'] == 'QuantumEspresso': TODO
        #    __init_qe_(self, filename)

        if general_par['verbosity'] > 0:
            if mpi.is_master_node():
                print("Element names:\n", self.struct.el_names)
                print("Q-types:\n", self.struct.q_types)
                print("Number of types:\n", self.struct.ntypes)
                print("Total number of atoms:\n", self.struct.nq)
                print("Number of each atom:\n", self.struct.nions)
                print("Type of ion:\n", self.struct.type_of_ion)
                print("Bravais Matrix:\n", self.struct.a_brav)

        counts = {}
        for i in self.struct.type_of_ion:
            counts[i] = counts.get(i, 0) + 1
        self.type_ion = [self.struct.el_names[i] + str(i) for i in counts for _ in range(counts[i])]
        # self.type_ion = [self.struct.el_names[i] for i in counts for _ in range(counts[i])]
        self.struct.q_types = np.concatenate(self.struct.q_types)

        if run_only == True:
            mpi.report("\nIon names for [MAGINT] :\n", self.type_ion)
            sys.exit()

        if self.single_atom:
            self.atom_list = [(self.atom0, count) for count in range(len(interact_sites[self.atom0]))]
        else:
            self.atom_list = [(atom, count) for atom in interacting_types for count in range(len(interact_sites[atom]))]

        for ind0 in self.atom_list:
            atom, count = ind0
            self.R_origin[ind0] = {}
            self.R_origin[ind0] = self.struct.q_types[self.type_ion.index(atom) + count]
            # for ind1 in range(len(self.R_origin[ind0])):
            #     if self.R_origin[ind0][ind1] > 0.5:
            #         self.R_origin[ind0][ind1] -= 1

        # Convert to Cartesian Coordinates
        # for i in range(len(self.struct.q_types)):
        #     self.struct.q_types[i, :] = self.struct.frac_to_carth(self.struct.q_types[i, :])
        # self.struct.q_types[i, :] = np.dot(self.struct.a_brav, self.struct.q_types[i, :].T).T

        self.latt_carth = self.struct.a_brav

        if general_par['dft_exec'] == 'Wien2k':
            filename = os.getcwd().rpartition('/')[2]
            f = filename + ".outputnn"
            self.nshells, self.R_shell, self.R_vec, self.T_vec, self.R_site, self.coor_num = self.__read_type_(
                f, max_n_shells, tol=tol_shells)
        else:
            self.nshells, self.R_shell, self.R_vec, self.T_vec, self.R_site, self.coor_num = self.find_nn(max_n_shells,
                                                                                                          direct=direct,
                                                                                                          tol=tol_shells)

        if general_par['verbosity'] > 0:
            mpi.report('-' * 40)
            mpi.report("\nCartesian coordinates of equivalent sites:")
            for atom in set(interacting_types):
                mpi.report("\nATOM: %s" % (atom))
            mpi.report("\nNum_shells = %s" % (self.nshells))
            for key0 in self.R_vec:
                print("\nIon Labels = ", key0)
                mpi.report("\nR_shells = %s" % (self.R_shell[key0]))
                for i in range(self.nshells):
                    if i >= len(self.coor_num[key0]):
                        mpi.report("\n\nThe range of Shells extends beyond the calculated nearest neighbors %s !\n")
                        break
                    mpi.report("\nShell = %s" % (i))
                    mpi.report("\nCoord_num = %s" % (self.coor_num[key0][i]))
                    if print_vecs == 'Basis' and general_par['verbosity'] > 1:
                        mpi.report("\nConnecting unit-cell vectors :\n")
                        for i1 in range(self.coor_num[key0][i]):
                            key = "%s_%s" % (i, i1)
                            mpi.report("%s :  %s" % (key, self.T_vec[key0][key]))
                        mpi.report("\nConnecting nearest-neighbors vectors :\n")
                        for i1 in range(self.coor_num[key0][i]):
                            key = "%s_%s" % (i, i1)
                            mpi.report("%s :  %s" % (key, self.R_vec[key0][key]))

    def __init_Vasp_(self, general_par):

        # Read File using dft_tools
        filename = 'POSCAR'
        struct = Poscar()
        path = general_par['folder'] + '/'
        struct.from_file(vasp_dir=path, poscar_filename=filename)
        self.struct = struct

    def __init_Wien2k_(self, filename):

        # Read File using parser
        filename = os.getcwd().rpartition('/')[2]
        f = filename + ".outputnn"
        struct = Wien2kStruct(f)
        self.struct = struct
        # self.struct.a_brav = self.struct.a_brav.T

    def calc_d_origin(self, verbosity=1):
        """
        Computes the minimal connecting vectors between the interaction sites
        of two atoms.

        Parameters
        ----------
        verbosity : int, optional
            Determines the amount of output information. Higher values produce more verbose output.
            Default value is 1.

        Returns
        -------
        d_origin : dict
            A dictionary storing minimal connecting vectors between interaction
            sites for the specified atom pairs. Keys are tuples in the form
            (atom_type1, site1, atom_type2, site2), and the corresponding values
            are the minimal connecting vectors. This dictionary is also set as an
            attribute of the object.

        Notes
        -----
        The function uses R_arr_set to determine if lattice vectors have been
        read from an output file. If not, R_origins are assumed to be [0,0,0].
        The function also uses the method __shift_vec_ to correct for vectors
        crossing periodic boundaries.
        """
        self.d_origin = {}
        atom0 = self.atom0
        atom1 = self.atom1
        if self.R_arr_set:
            # if lattice vectors are not read from outputnn, R_origins are [0,0,0]
            for i in self.interact_sites[self.atom0]:
                for j in self.interact_sites[self.atom1]:
                    self.d_origin[((atom0, i), (atom1, j))] = array([0.0, 0.0, 0.0])
        else:
            if verbosity > 0:
                mpi.report("\n\nMinimal connecting vectors between sites of types %s and %s:" % (atom0, atom1))
            for ind0 in self.atom_list:
                if not self.calc_all_sites and ind0 != self.atom_list[0]: continue
                # ind0 = self.atom_list[0]
                atom0, count0 = ind0
                for ind1 in self.atom_list:
                    atom1, count1 = ind1
                    dd = self.R_origin[ind0] - self.R_origin[ind1]
                    self.d_origin[((atom0, count0), (atom1, count1))] = dd
                    if verbosity > 0:
                        mpi.report("\n%s site-%s -> %s site-%s:  %7.3f  %7.3f %7.3f" % (
                            atom0, count0, atom1, count1, dd[0], dd[1], dd[2]))

    def __shift_vec_(self, dd, ir=1):
        """
        tries to reduce the distance between sites by applying translations

        Inputs
        -------
        Interacting types :
        Interacting sites :
        Interacting max shells :
        Tol :

        Returns
        -------

        """
        len_in = abs(dot(dd, dd))
        for ix in range(-ir, ir + 1, 1):
            for iy in range(-ir, ir + 1, 1):
                for iz in range(-ir, ir + 1, 1):
                    lvec = dot(array([ix, iy, iz]), self.latt_carth)
                    dd_new = dd + lvec
                    len_new = abs(dot(dd_new, dd_new))
                    if (len_new < len_in):
                        dd = dd_new
                        len_in = len_new
        return dd

    def find_nn(self, max_n_shells, direct, tol):
        """
        Finds the nearest neighbors, next nearest neighbors, and so on from the POSCAR
        file in VASP or case.outputnn in Wien2k. It also initializes the structural data necessary for MagInt.

        Parameters
        ----------
        max_n_shells : int
            Maximum number of shells to be considered when finding neighbors.
        direct : bool
            If True, returns the basis vector. If False, returns the dot product
            of the a_brav.T matrix and the basis vector.

        Returns
        -------
        nshells : int
            Number of shells found.
        R_shell : list of floats
            Sorted list of unique distances, indicating the radial distances
            for each shell.
        R_vec : dict
            Dictionary mapping each shell and coordinate number to the corresponding
            vector. The keys are strings in the format "nshell_coorNum".
        R_site : dict
            Dictionary similar to R_vec, but the values are initialized to 0.
        coor_num : list of ints
            A list indicating the coordination number for each shell.

        Notes
        -----
        The function uses the atom_list attribute to determine if only one atom type
        is present. It employs a brute-force approach, computing distances for all
        possible shell combinations up to max_n_shells, and then sorts and filters
        the results to extract the relevant shells and coordination numbers.
        """
        # Initialize
        R_vec = {}
        T_vec = {}
        R_site = {}
        R_shells = {}
        C_num = {}
        # diff_vec = []
        tol = int(-np.floor(np.log10(tol)))

        shift = list(itertools.product(range(-max_n_shells - 1, max_n_shells + 1), repeat=3))

        for ind0 in self.atom_list:
            if not self.calc_all_sites and ind0 != self.atom_list[0]: continue
            for ind1 in self.atom_list:

                R_bas_vec = []
                distance = []
                coor_num = []
                R_vec[(ind0, ind1)] = {}
                T_vec[(ind0, ind1)] = {}
                R_site[(ind0, ind1)] = {}
                R_shells[(ind0, ind1)] = {}
                nshells = 0

                diff_vec = self.R_origin.get(ind0) - self.R_origin.get(ind1)

                for i in range(len(shift)):
                    # for j in range(len(diff_vec)):
                    # T_vec_ini = np.dot(self.struct.a_brav, np.array(shift[i] + diff_vec).T).T  # [j]
                    T_vec_ini = np.dot(self.struct.a_brav, np.array(shift[i]).T).T + diff_vec  # [j]
                    dist = np.round(np.linalg.norm(T_vec_ini), tol)
                    if dist > 0:
                        distance.append(dist)
                        R_bas_vec.append(shift[i])

                sort_bas = sorted(zip(distance, R_bas_vec), key=lambda x: x[0])
                distance = sorted(distance)
                R_shell = sorted(list(set([x[0] for x in sort_bas])))
                R_bas_vec = [x[1] for x in sort_bas]

                count = 0
                coor_num.append(0)

                for count1, i in enumerate(distance):
                    if i == R_shell[count]:
                        coor_num[count] += 1
                    else:
                        coor_num.append(0)
                        nshells += 1
                        count += 1
                        coor_num[count] += 1
                    if nshells >= max_n_shells:
                        break
                    str1 = "%s_%s" % (nshells, coor_num[count] - 1)
                    R_site[(ind0, ind1)][str1] = 0
                    if direct:
                        T_vec[(ind0, ind1)][str1] = np.array(R_bas_vec[count1])
                        # R_vec[str1] = R_bas_vec[count1] + diff_vec
                        R_vec[(ind0, ind1)][str1] = np.round((self.struct.a_brav @ R_bas_vec[count1]) + diff_vec, 7)
                    else:
                        T_vec[(ind0, ind1)][str1] = np.round(self.struct.a_brav @ R_bas_vec[count1], 7)
                        R_vec[(ind0, ind1)][str1] = np.round((self.struct.a_brav @ R_bas_vec[count1]) + diff_vec, 7)

                R_shell = R_shell[:nshells]
                R_shells[(ind0, ind1)] = R_shell
                C_num[(ind0, ind1)] = coor_num

        return nshells, R_shells, R_vec, T_vec, R_site, C_num

    def __read_type_(self, f, max_n_shells, tol=2e-4):
        """
        Reads and processes atomic neighbor data from a Wien2k case.outputnn file. This function analyzes
        the structural data of specified atom types, calculates the nearest neighbor vectors, and categorizes
        them into shells based on their distances.

        Parameters
        ----------
        f : str or file object
            Filename or file object from which to read the atomic neighbor data.
        max_n_shells : int
            Maximum number of neighbor shells to consider when categorizing the atomic vectors.
        tol : float, optional
            Tolerance value to distinguish between different shells. Atoms within this distance from each other
            are considered in the same shell. Default is 2e-4.

        Returns
        -------
        max_n_shells : int
            The number of shells actually found, which may be less than or equal to the input 'max_n_shells'.
        R_shells : dict
            Dictionary where keys are tuples of atom pairs and values are lists of unique shell distances for
            those atom pairs.
        R_vec : dict
            Dictionary mapping atom pairs to their relative position vectors, categorized by shells.
        T_vec : dict
            Dictionary mapping atom pairs to vectors relative to the translation vectors without basis vectors,
            categorized by shells.
        R_site : dict
            Dictionary that maps shells and coordination numbers to equivalent sites for atom pairs.
        C_num : dict
            Dictionary that keeps track of coordination numbers for each atom pair in different shells.

        Raises
        ------
        ValueError
            If the file does not contain the expected atomic data or format.

        """
        # Read lattice parameters
        R_vec = {}
        T_vec = {}
        R_site = {}
        R_shells = {}
        C_num = {}

        f = open(f)

        for ind0 in self.atom_list:
            if not self.calc_all_sites and ind0 != self.atom_list[0]: continue
            # ind0 = self.atom_list[0]
            for ind1 in self.atom_list:
                if not self.calc_all_sites and ind1 == ind0 and len(self.atom_list) > 1: continue
                # Retrieve the actual name, position of each element
                atom0 = self.struct.el_names[self.type_ion.index(ind0[0])]
                num0 = int(self.type_ion.index(ind0[0]))
                atom1 = self.struct.el_names[self.type_ion.index(ind1[0])]
                num1 = int(self.type_ion.index(ind1[0]))
                # atom0 = ''.join(filter(str.isalpha, ind0[0]))

                coor_num = []
                R_vec[(ind0, ind1)] = {}
                T_vec[(ind0, ind1)] = {}
                R_site[(ind0, ind1)] = {}
                R_shells[(ind0, ind1)] = {}
                R_shell = []
                nshells = -1
                f.seek(0)

                for il in f:
                    if (il.find('ATOM:') == 1 and il.find(atom0) > -1 and il.find('EQUIV.') > -1 and
                            int(il.split()[il.split().index('EQUIV.') + 1]) - 1 == ind0[1] and int(
                                il.split()[il.split().index('ATOM:') + 1]) - 1 == num0):
                        for ist in range(len(il.split())):
                            if il.split()[ist] == 'AT': break
                        coor = array([float(il.split()[i]) for i in range(ist + 1, ist + 4)])
                        dd = 0.0
                        # include lattice vectors only for the inequivalent site for each corr. type
                        # if (n_eq_site0 != 0): continue
                        for il1 in f:
                            if (il1.find('ATOM:') == 1 and il1.split()[2] == atom1 and il1.find('ANG') > -1 and int(
                                    il1.split()[-1]) - 1 == ind1[1] and int(
                                il1.split()[il.split().index('ATOM:') + 1]) - 1 == num1):
                                for ist in range(len(il1.split())):
                                    if il1.split()[ist] == 'AT': break
                                #coor1 = array([float(il1.split()[i]) for i in range(ist + 1, ist + 4)])
                                x=float(il1[23:31])
                                y=float(il1[31:39])
                                z=float(il1[39:47])
                                coor1 = array([x,y,z])
                                #coor1 = array([float(il1.split()[i]) for i in range(ist + 1, ist + 4)])
                                #n_eq_site1 = int(il1.split()[ist + 9]) - 1
                                n_eq_site1 = int(il1[80:85])
                                coor1 = coor1 - coor

                                diff_vec = self.R_origin.get(ind1) - self.R_origin.get(ind0)

                                dd_new = sqrt(dot(self.struct.frac_to_carth(coor1), self.struct.frac_to_carth(coor1)))

                                if (abs(dd - dd_new) > tol):
                                    if (nshells == max_n_shells - 1):
                                        nshells += 1
                                        R_shells[(ind0, ind1)] = R_shell[:nshells]
                                        C_num[(ind0, ind1)] = coor_num
                                        break
                                        # return nshells, R_shells, R_vec, T_vec, R_site, C_num
                                    nshells += 1
                                    R_shell.append(dd_new)
                                    coor_num.append(1)
                                    dd = dd_new
                                else:
                                    coor_num[-1] += 1
                                str1 = "%s_%s" % (nshells, coor_num[-1] - 1)
                                R_site[(ind0, ind1)][str1] = n_eq_site1
                                T_vec[(ind0, ind1)][str1] = self.struct.frac_to_carth(np.round(coor1 - diff_vec, 3))
                                R_vec[(ind0, ind1)][str1] = self.struct.frac_to_carth(coor1)
                            elif (il1.find('ATOM:') == -1 and il1.find('RMT(') == -1 and il1.find('SUMS') == -1):
                                break

                if R_shell:
                    R_shells[(ind0, ind1)] = R_shell[:max_n_shells]
                C_num[(ind0, ind1)] = coor_num

                nshells += 1

        return max_n_shells, R_shells, R_vec, T_vec, R_site, C_num
