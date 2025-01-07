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
from MagInt.utils import *
from math import cos, sin, sqrt, radians
import numpy as np
import sys


class Wien2kStruct:
    """
    The Wien2kStruct class initializes the structure file as in dft_tools

    """

    def __init__(self, filename):
        self.filename = filename
        self.atom_positions = {}
        self.q_types = []
        self.el_names = []

        # New attributes
        self.ntypes = 0
        self.nq = 0
        self.nions = []
        self.type_of_ion = []

        # Attributes for lattice vectors and angles
        self.lattice_vecs = None
        self.angles = None
        # Bravais matrix
        self.bravais = None
        self.a_brav = None
        self.orth = True
        self.lattype = None

        self.parse_struct_file()

    def parse_struct_file(self):
        """
        Parses the input file to extract lattice information, atom positions, and other
        relevant structural data. It updates the class attributes accordingly.

        Attributes Modified
        -------------------
        lattype : str
            Lattice type (e.g., "P", "C", etc.).
        atom_positions : dict
            Dictionary mapping atom names to their positions.
        el_names : list
            List of element names found in the file.
        q_types : list
            List of atomic positions for each element type.
        nions : list
            List of the number of ions for each element type.
        lattice_vecs : ndarray
            Lattice vectors.
        angles : ndarray
            Lattice angles.
        bravais : ndarray
            The Bravais matrix for the structure.
        ntypes : int
            The number of unique atom types.
        nq : int
            Total number of atoms.
        type_of_ion : list
            Type of ion corresponding to each atom in the system.
        orth : bool
            Whether the system is orthorhombic.
        """
        with open(self.filename, 'r') as file:
            lines = file.readlines()

            self.lattype = lines[1].strip()[0:]

            current_atom_id = None
            current_positions = []
            found_rela = False

            for i, line in enumerate(lines):
                split_line = line.strip().split()
                cs1 = split_line and split_line[0].startswith('-') and split_line[0][1:].isdigit()
                cs2 = split_line and len(split_line) == 4 and split_line[0].isdigit()
                # if split_line and split_line[0].startswith('-') and split_line[0][1:].isdigit():
                if cs1 or cs2:
                    identifier = int(split_line[0])

                    if current_atom_id is not None and identifier != current_atom_id:
                        current_positions = []

                    current_atom_id = identifier
                    coords = [float(val) for val in split_line[1:]]
                    current_positions.append(coords)

                elif "NPT" in line:
                    atom_name = line.split()[0]
                    self.atom_positions[atom_name] = current_positions
                    self.el_names.append(atom_name)
                    self.q_types.append(np.array(current_positions))
                    self.nions.append(len(current_positions))
                    current_positions = []

                elif "RELA" in line:
                    found_rela = True
                    width = 10
                    rela_line = lines[i + 1].strip()
                    self.lattice_vecs = np.array(
                        [float(rela_line[j:j + width].strip()) for j in range(0, 3 * width, width)])
                    # commented line does not work for 3-digit angles
                    # self.angles = np.array([float(x) for x in rela_data[3:]])
                    ang = []
                    for p in [30, 40, 50]: ang.append(float(lines[i + 1][p:p + 10]))
                    self.angles = np.array(ang)

                elif "Bravais Matrix:" in line:
                    bravais_matrix = []
                    for j in range(3):
                        bravais_row = [float(x) for x in lines[i + j + 1].strip().split()]
                        bravais_matrix.append(bravais_row)
                    self.bravais = np.array(bravais_matrix)

                elif "Bond-Valence" in line:
                    break

            # Calculate the new attributes
            self.ntypes = len(self.el_names)
            self.nq = sum(self.nions)
            for i, count in enumerate(self.nions):
                self.type_of_ion.extend([i] * count)

        for i in range(3):
            if abs(self.angles[i] - 90.0) > 1e-6: self.orth = False
        for i in range(3): self.angles[i] = radians(self.angles[i])

        self.calc_a_brav()

    def calc_a_brav(self):
        """
        Calculates the Bravais lattice vectors in Cartesian coordinates
        based on the provided Bravais matrix and lattice type. The resulting
        lattice vectors are stored in the `a_brav` attribute of the class.

        Attributes Modified
        -------------------
        a_brav : ndarray
            The Bravais lattice vectors in Cartesian coordinates.
        """
        self.a_brav = zeros((3, 3))

        # first find lattice translations in carthesian coords.
        if self.lattype[0:1] != 'R':
            # self.a_brav = (M @ self.bravais).T
            self.a_brav[0, :] = self.frac_to_carth(array([1.0, 0.0, 0.0]))
            self.a_brav[1, :] = self.frac_to_carth(array([0.0, 1.0, 0.0]))
            self.a_brav[2, :] = self.frac_to_carth(array([0.0, 0.0, 1.0]))
        else:
            # for R lattices use hexagonal cell
            self.a_brav[0, :] = array([0.8660254 * self.lattice_vecs[0], -0.5 * self.lattice_vecs[1], 0.0])
            self.a_brav[1, :] = array([0.0, self.lattice_vecs[1], 0.0])
            self.a_brav[2, :] = array([0.0, 0.0, self.lattice_vecs[2]])
        return

    def frac_to_carth(self, coor):
        """
        Converts fractional coordinates to Cartesian coordinates based on
        the given lattice vectors and Bravais matrix. The conversion considers
        different lattice types and angles.

        Parameters
        ----------
        coor : ndarray
            Fractional coordinates to be converted.

        Returns
        -------
        ndarray
            Cartesian coordinates corresponding to the given fractional coordinates.
        """
        if self.orth == True:
            # a,b,c orthogonal
            coor_car = coor * self.lattice_vecs
        else:
            if self.lattype[0:1] == 'R' or self.lattype[0:1] == 'H':
                coor_car = np.dot(self.bravais, coor) * self.lattice_vecs
                # coor_car = coor * self.lattice_vecs
            elif self.lattype[0:3] == 'CXZ':
                singam = sin(self.angles[2])
                cosgam = cos(self.angles[2])
                coor_car = zeros((3,))
                coor_car[0] = coor[0] * singam * self.lattice_vecs[0]
                coor_car[1] = coor[0] * cosgam * self.lattice_vecs[0] + coor[1] * self.lattice_vecs[1]
                coor_car[2] = coor[2] * self.lattice_vecs[2]
            else:
                # coor_car = coor * self.lattice_vecs
                coor_car = dot(self.bravais, coor * self.lattice_vecs)
                # coor_car = dot(self.bravais.transpose(), coor * self.lattice_vecs)
                # coor_car = self.latt_carth.transpose() @ coor
                # coor_car=dot(coor*self.lpar,self.br2)
        return coor_car
