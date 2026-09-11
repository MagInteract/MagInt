.. highlight:: bash
.. _installation:

Installation 
#############

Prerequisites
-------------

#. The `TRIQS` library version 4.0.x, see `TRIQS installation instruction <https://triqs.github.io/triqs/latest/install.html>`_.

#. `TRIQS/DFTTools <https://triqs.github.io/dft_tools/latest/>`_, together with
   `TRIQS/dftkit <https://triqs.github.io/dftkit/latest/>`_, which provides the
   DFT converters from TRIQS 4.0 onwards. Both are imported at run time.

#. A Fortran compiler and LAPACK, used to build the Hubbard-I solver.

#. Make sure to install besides the triqs requirements also the python packages::

     $ pip3 install --user scipy pytest f90wrap meson ninja

   ``f90wrap`` generates the Fortran wrapper at configure time and is imported
   again at run time by the generated module, so it has to remain installed.
   ``meson`` and ``ninja`` are the build backend that ``f2py`` uses from NumPy
   1.26 onwards.

#. To build the documentation the following extra python packages are needed::

     $ pip3 install --user sphinx sphinx-autobuild pandoc nbsphinx linkify-it-py sphinx_rtd_theme myst-parser



Manual installation via CMake
-----------------------------

We provide hereafter the build instructions in the form of a documented bash script. Please change the variable
INSTALL_PREFIX to point to your TRIQS installation directory::
    
    INSTALL_PREFIX=/path/to/triqs
    # source the triqsvars.sh file from your TRIQS installation to load the TRIQS environment
    source $INSTALL_PREFIX/share/triqs/triqsvars.sh

    # clone the MagInt repository from GitHub
    git clone https://github.com/MagInteract/MagInt.git magint.src

    # Create and move to a new directory where you will compile the code
    mkdir magint.build && cd magint.build

    # In the build directory call cmake, including any additional custom CMake options, see below
    cmake ../magint.src

    # Compile the code, run the tests, and install the application
    make
    make test
    make install

This installs MagInt into your TRIQS installation folder.

To build ``MagInt`` with documentation you should run::

     $ cmake path/to/magint.src -DBUILD_DOC=ON
     $ make 
     $ sphinx-autobuild path/to/magint.src/doc ./doc/html -c ./doc/


Version compatibility
---------------------

The release version ``MagInt`` 3.0 is compatible with TRIQS 4.0.x

Custom CMake options
--------------------

The compilation of ``MagInt`` can be configured using CMake-options::

    cmake ../magint.src -DOPTION1=value1 -DOPTION2=value2 ...

+-----------------------------------------------------------------+-----------------------------------------------+
| Options                                                         | Syntax                                        |
+=================================================================+===============================================+
| Specify an installation path other than path_to_triqs           | -DCMAKE_INSTALL_PREFIX=path_to_magint         |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the documentation                                         | -DBUILD_DOC=ON                                |
+-----------------------------------------------------------------+-----------------------------------------------+

``CMAKE_BUILD_TYPE`` is not configurable on the command line: ``CMakeLists.txt``
sets it to ``Release`` unconditionally, which overrides any value passed to
``cmake``. The test suite is likewise always configured, so there is no option
to disable it.
