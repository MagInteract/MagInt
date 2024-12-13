![logo_MagInt](https://github.com/dariofiosca/MagInt/blob/main/logo.png)


# __MagInt__

---

Copyright (C) 2012-2024: Leonid V. Pourovskii, 2023-2024: Dario Fiore Mosca

This program implements the force-theorem Hubbard-I (FT-HI) approach to calculation of intersite exchange interactions in correlated insulators. The formalism is given in L. V. Pourovskii Phys. Rev. B 94, 115117 (2016). It is a post-processing program that works on top of a fully converged DFT+HI calculation, and for this purpose it is interfaced with the [TRIQS](https://triqs.github.io/triqs/latest/) software library, and the DFT code interface [TRIQS/DFTTools](https://triqs.github.io/dft_tools/latest/). 

MagInt has been interfaced with the following bandstructure codes: [Wien2k](http://www.wien2k.at) version 13 to 15.  [VASP](https://www.vasp.at) version 6.3.0 or above.

MagInt works with python < 3.12. 

### Documentation & tutorials

To learn how to use MagInt [Read the Manual](docs/magint_manual.pdf)


### Installation

You can install MagInt via CMake

```
rm -rf MagInt.build

# Set this variable to your desired install directory
INSTALL_PREFIX=$(pwd)/install

# Set the number of cores for the compilation
NCORES=4

# Use cmake to configure the maginteract build process
mkdir -p MagInt.build && cd MagInt.build

cmake ../MagInt.src -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX

# Build, test and install triqs
make -j $NCORES
make install
make test
cd ../
```

### Usage

In the following the main commands are outlined. For each step a 'maginteract.ini' file is required. Have a look at the Documentation and Tutorials for more 
information 

1. Calculation Density of States

```
python /path/to/MagInt/main.py dos
```

2. Shifting the chemical potential

```
python /path/to/MagInt/main.py shift_mu
```

3. Calcualting the Pseudo-Spin basis

For this step a problem-specific script is required.

4.  Calculation of Intersite-exchange-interactions

```
python /path/to/MagInt/main.py magint 
```

---

Please make sure that you have a valid TRIQS and TRIQS/DFTTools installation matching the required versions for MagInt. 

Within the current version of MagInt we provide FORTRAN Hubbard-I solver. 

---

authors: L. V. Pourovskii and D. Fiore Mosca

### LICENCE

This application is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version (see http://www.gnu.org/licenses/).

It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
