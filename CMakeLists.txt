# ------------------------------------------------------------------------------------#
# MagInt library
# Written by  : Leonid V. Pourovskii (CPHT Ecole Polytechnique) 2012-2024
#             : Dario Fiore Mosca (CPHT Ecole Polytechnique) 2023-2024
# Email: leonid@cpht.polytechnique.fr
# ------------------------------------------------------------------------------------#
#
#    MagInteract implements the force-theorem Hubbard-I (FT-HI) approach to intersite
#    exchange interactions in correlated insulators.
#    The formalism is given in L. V. Pourovskii Phys. Rev. B 94, 115117 (2016)
#    This Python-3 version is based on TRIQS library
#
# ------------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
cmake_policy(VERSION 3.20)
if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif()

# ############
# Define Project
project(MagInt VERSION 3.0.0 LANGUAGES CXX Fortran)
get_directory_property(IS_SUBPROJECT PARENT_DIRECTORY)

# Version number of the application
set (MAGINT_VERSION "3.0")
set (MAGINT_RELEASE "3.0.0")

# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/cmake/Modules)

# start configuration
set(CMAKE_BUILD_TYPE Release)
enable_testing()

# We use shared libraries
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)

# Load PYTHON
find_package(Python3 COMPONENTS Interpreter Development)

if (NOT ${TRIQS_WITH_PYTHON_SUPPORT})
 MESSAGE(FATAL_ERROR "Hubbard I application require Python support in TRIQS")
endif()

# Check that versions are compatible
if(NOT MAGINT_VERSION EQUAL TRIQS_VERSION)
  message(FATAL_ERROR "The application version is not compatible with the TRIQS library (TRIQS library version: ${TRIQS_VERSION} while this application version: ${MAGINT_VERSION})")
endif()


# Default Install directory to TRIQS_ROOT if not given or invalid.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS "No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
  set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT FALSE)
endif()
if(NOT IS_SUBPROJECT)
  message(STATUS "-------- CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} --------")
endif()
set(${PROJECT_NAME}_BINARY_DIR ${PROJECT_BINARY_DIR} CACHE STRING "Binary directory of the ${PROJECT_NAME} Project")

# Codes and tests
add_subdirectory(fortran/MagInt)
add_subdirectory(python/MagInt)
add_subdirectory(test)

# Documentation
option(BUILD_DOC "Build documentation" OFF)
if(${BUILD_DOC})
  if(NOT TRIQS_WITH_DOCUMENTATION)
    message("Error: TRIQS library has not been compiled with its documentation")
  endif()
  add_subdirectory(doc)
endif(${BUILD_DOC})



