.. highlight:: bash
.. _installation:

Installation 
#############

Prerequisites
-------------

#. The :ref:`TRIQS <triqslibs:welcome>` library, see `TRIQS installation instruction <https://triqs.github.io/triqs/latest/install.html>`_.

#. Make sure to install besides the triqs requirements also the python packages::

     $ pip3 install --user scipy pytest f90wrap

#. To build the documentation the following extra python packages are needed::

     $ pip3 install --user sphinx sphinx-autobuild pandoc nbsphinx linkify-it-py sphinx_rtd_theme myst-parser



Manual installation via CMake
-----------------------------

We provide hereafter the build instructions in the form of a documented bash script. Please change the variable
INSTALL_PREFIX to point to your TRIQS installation directory::
    
    INSTALL_PREFIX=/path/to/triqs
    # source the triqsvars.sh file from your TRIQS installation to load the TRIQS environment
    source $(INSTALL_PREFIX)/share/triqs/triqsvars.sh

    # clone the MagInt repository from GitHub
    git clone https://github.com/MagInteract/MagInt.git magint.src

    # Create and move to a new directory where you will compile the code
    mkdir magint.build && cd magint.build

    # In the build directory call cmake, including any additional custom CMake options, see below
    cmake ../magint.src

    # Compile the code, run the tests, and install the application
    make test
    make install

This installs MagInt into your TRIQS installation folder.

To build ``MagInt`` with documentation you should run::

     $ cmake path/to/magint.src -DBuild_Documentation=ON
     $ make 
     $ sphinx-autobuild path/to/magint.src/doc ./doc/html -c ./doc/


Version compatibility
---------------------

The release version ``MagInt`` 1.0 is compatible with TRIQS 3.1.x

Custom CMake options
--------------------

The compilation of ``MagInt`` can be configured using CMake-options::

    cmake ../magint.src -DOPTION1=value1 -DOPTION2=value2 ...

+-----------------------------------------------------------------+-----------------------------------------------+
| Options                                                         | Syntax                                        |
+=================================================================+===============================================+
| Specify an installation path other than path_to_triqs           | -DCMAKE_INSTALL_PREFIX=path_to_solid_dmft     |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build in Debugging Mode                                         | -DCMAKE_BUILD_TYPE=Debug                      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Disable testing (not recommended)                               | -DBuild_Tests=OFF                             |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the documentation                                         | -DBuild_Documentation=ON                      |
+-----------------------------------------------------------------+-----------------------------------------------+