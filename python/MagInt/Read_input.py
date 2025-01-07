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
import configparser
import ast

# Define permitted configuration
general_keys = ['dft_exec', 'name', 'verbosity', 'folder', 'dft_filename', 'phase', 'tol_shells', 'use_symmetries']
solver_keys = ['beta', 'u_int', 'j_hund', 'eal_read', 'use_spin_orbit',
               'n_lev', 'n_moments', 'ommin', 'ommax', 'n_omega', 'dc_type',
               'kanamori', 'chemical_potential_init', 'ls', 'n_iomega', 'n_el', 'broadening', 'tol_eal']
basis_keys = ['split', 'ps_j', 'action']
magint_keys = ['atom1', 'atom2', 'n_shells', 'correlated_atoms', 'mult_conv', 'calc_all_sites','included_atoms']

# Define default configuration
default_config = configparser.ConfigParser()
default_config['GENERAL'] = {'dft_exec': "'Wien2k'", 'verbosity': 1, 'folder': '.', 'dft_filename': 'test',
                             'phase': 1, 'tol_shells': 0.001, 'use_symmetries': True}
default_config['SOLVER'] = {'beta': 40.0, 'u_int': 0.0, 'j_hund': 0.0, 'eal_read': False, 'use_spin_orbit': True,
                            'n_lev': 4, 'n_moments': 4, 'ommin': 0.0, 'n_omega': 2050, 'ls': 2,
                            'dc_type': 0, 'kanamori': True, 'chemical_potential_init': 0.0,
                            'n_iomega': 1025, 'n_el': 1, 'tol_eal': 1e-5}
default_config['BASIS'] = {'split': 0.000001, 'ps_j': 0, 'action': "'merge'"}
default_config['MAGINT'] = {'atom1': 'Atom1', 'atom2': 'Atom2', 'n_shells': 3, 'mult_conv': "'sph_tensor'",
                            'calc_all_sites': False, 'included_atoms': 0}


def read_input_file(filename):
    config = configparser.ConfigParser()
    config.read(filename)
    general_par = {'dft_exec': str, 'verbosity': int, 'name': str, 'folder': str, 'dft_filename': str, 'phase': float,
                   'tol_shells': float, 'use_symmetries': bool}
    solver_par = {'beta': float, 'u_int': list, 'j_hund': float, 'eal_read': bool, 'use_spin_orbit': bool, 'n_lev': int,
                  'n_moments': int, 'ommin': float, 'ommax': float, 'n_omega': int, 'dc_type': int,
                  'kanamori': bool, 'chemical_potential_init': list, 'ls': list, 'n_iomega': int, 'n_el': list,
                  'broadening': float, 'tol_eal': float}
    basis_par = {'split': float, 'ps_j': list, 'action': str}
    magint_par = {'atom1': str, 'atom2': str, 'n_shells': int, 'correlated_atoms': list, 'mult_conv': str,
                  'calc_all_sites': bool, 'included_atoms': list}

    for section in default_config.keys():
        if config.has_section(section):
            for key in default_config[section].keys():
                if key not in config[section]:
                    config[section][key] = default_config[section][key]

    mpi.report('-' * 40)
    mpi.report('-' * 40 + '\nGeneral Parameters:\n')
    for key, value in config['GENERAL'].items():
        if key != 'phase':  # DEBUG
            mpi.report('{0: <20} {1: <4}'.format(key, str(value)))
        if key not in general_par:
            raise ValueError(f'{key} is not a valid general parameter name')
        if key == 'name' and value == 'default':
            value = os.getcwd().rpartition('/')[2]
        general_par[key] = ast.literal_eval(value)

    if 'SOLVER' in config.sections():
        mpi.report('-' * 40 + '\nSolver Parameters:\n')
        for key, value in config['SOLVER'].items():
            mpi.report('{0: <20} {1: <4}'.format(key, str(value)))
            if key not in solver_par:
                raise ValueError(f'{key} is not a valid solver parameter name')
            solver_par[key] = ast.literal_eval(value)

    if 'BASIS' in config:
        mpi.report('-' * 40 + '\nBasis Parameters\n')
        for key, value in config['BASIS'].items():
            mpi.report('{0: <20} {1: <4}'.format(key, str(value)))
            if key not in basis_par:
                raise ValueError(f'{key} is not a valid basis parameter name')
            basis_par[key] = ast.literal_eval(value)

    if 'MAGINT' in config:
        mpi.report('-' * 40 + '\nMagInt Parameters\n')
        for key, value in config['MAGINT'].items():
            mpi.report('{0: <20} {1: <4}'.format(key, str(value)))
            if key not in magint_par:
                raise ValueError(f'{key} is not a valid basis parameter name')
            magint_par[key] = ast.literal_eval(value)

    return general_par, solver_par, basis_par, magint_par


