#triqs_build_f2py_module( LP_hubbardI hubbard_I hubbard_I.pyf hubbard_I.f90)

find_package(LAPACK)
if (NOT LAPACK_FOUND)
    message(FATAL_ERROR "LAPACK not found!")
endif()

# f90wrap ships its drivers as Python console scripts. A `pip install --user`
# puts those in a script directory that is usually not on PATH, so ask the
# interpreter we were configured with where they would live and search there too.
# Locating them up front makes a missing installation fail here with an
# actionable message, rather than letting the steps below silently produce nothing.
execute_process(COMMAND ${Python3_EXECUTABLE} -c
                  "import sysconfig; print(sysconfig.get_path('scripts','posix_user'), sysconfig.get_path('scripts'), sep=';')"
                OUTPUT_VARIABLE PYTHON_SCRIPT_DIRS
                OUTPUT_STRIP_TRAILING_WHITESPACE
                COMMAND_ERROR_IS_FATAL ANY)

find_program(F90WRAP_EXECUTABLE f90wrap HINTS ${PYTHON_SCRIPT_DIRS})
find_program(F2PY_F90WRAP_EXECUTABLE f2py-f90wrap HINTS ${PYTHON_SCRIPT_DIRS})
if(NOT F90WRAP_EXECUTABLE OR NOT F2PY_F90WRAP_EXECUTABLE)
  message(FATAL_ERROR
    "f90wrap/f2py-f90wrap not found.\n"
    "Searched PATH and: ${PYTHON_SCRIPT_DIRS}\n"
    "Install with:  ${Python3_EXECUTABLE} -m pip install f90wrap")
endif()
message(STATUS "Found f90wrap: ${F90WRAP_EXECUTABLE}")

# FindLAPACK reports full paths to library files, which cannot be handed to
# f2py as-is. Split the first entry into a -L directory and a -l name.
list(GET LAPACK_LIBRARIES 0 _lapack_lib)
get_filename_component(_lapack_dir ${_lapack_lib} DIRECTORY)
get_filename_component(_lapack_name ${_lapack_lib} NAME_WE)
string(REGEX REPLACE "^lib" "" _lapack_name "${_lapack_name}")

execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -c -g -O3 -fPIC ${CMAKE_CURRENT_SOURCE_DIR}/hubbard_I.f90
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMAND_ERROR_IS_FATAL ANY)

execute_process(COMMAND ${F90WRAP_EXECUTABLE} -v -m hubbard_I ${CMAKE_CURRENT_SOURCE_DIR}/hubbard_I.f90 -k ${CMAKE_CURRENT_SOURCE_DIR}/kind_map
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMAND_ERROR_IS_FATAL ANY)

# numpy's meson f2py backend ignores -I when searching for Fortran .mod files,
# so the module directory has to reach the compiler through FFLAGS as well.
execute_process(COMMAND ${CMAKE_COMMAND} -E env FFLAGS=-I${CMAKE_CURRENT_BINARY_DIR}
                  ${F2PY_F90WRAP_EXECUTABLE} -c
                  -L${_lapack_dir} -l${_lapack_name}
                  -I${CMAKE_CURRENT_BINARY_DIR}
                  -m _hubbard_I
                  ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_hubbard_I.f90
                  ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_toplevel.f90
                  ${CMAKE_CURRENT_BINARY_DIR}/hubbard_I.o
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMAND_ERROR_IS_FATAL ANY)

execute_process(COMMAND python3-config --extension-suffix OUTPUT_VARIABLE SUFFIX OUTPUT_STRIP_TRAILING_WHITESPACE)

string(CONCAT FILE_NAME "_hubbard_I" ${SUFFIX})
# Export to the parent scope so python/MagInt can link against the built module.
set(HUBBARD_I_MODULE ${FILE_NAME} CACHE INTERNAL "name of the compiled _hubbard_I extension")

file(READ "${CMAKE_CURRENT_BINARY_DIR}/hubbard_I.py" file_content)
# Replace the line
string(REPLACE "import _hubbard_I" "from . import _hubbard_I" modified_content "${file_content}")
# Write the modified content back to the file
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/hubbard_I.py" "${modified_content}")

# The python sources import the generated wrapper as `hubbardI.hubbard_I`, so
# install it and its extension module as a `hubbardI` package. This mirrors the
# layout the build tree already exposes through the symlink made in python/MagInt.
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/__init__.py "")
set(PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/hubbardI)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/${FILE_NAME} DESTINATION ${PYTHON_LIB_DEST})
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/hubbard_I.py DESTINATION ${PYTHON_LIB_DEST})
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/__init__.py DESTINATION ${PYTHON_LIB_DEST})
