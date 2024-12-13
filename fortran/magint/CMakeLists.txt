#triqs_build_f2py_module( LP_hubbardI hubbard_I hubbard_I.pyf hubbard_I.f90)

execute_process(COMMAND gfortran -c -g -O3 -fPIC ${CMAKE_CURRENT_SOURCE_DIR}/hubbard_I.f90 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
execute_process(COMMAND f90wrap -v -m hubbard_I ${CMAKE_CURRENT_SOURCE_DIR}/hubbard_I.f90 -k ${CMAKE_CURRENT_SOURCE_DIR}/kind_map WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
find_package(LAPACK)
if (NOT LAPACK_FOUND)
    message(FATAL_ERROR "LAPACK not found!")
endif()

#execute_process(COMMAND sed -i '2s/import _hubbard_I/from . import _hubbard_I/' ${CMAKE_CURRENT_BINARY_DIR}/hubbard_I.py)
#message("Modifying Python import statement...")

execute_process(COMMAND f2py-f90wrap -c -L${LAPACK_LIBRARIES} -I${CMAKE_CURRENT_BINARY_DIR} -llapack -m _hubbard_I ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_hubbard_I.f90 ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_toplevel.f90 ${CMAKE_CURRENT_BINARY_DIR}/hubbard_I.o WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
execute_process(COMMAND python3-config --extension-suffix OUTPUT_VARIABLE SUFFIX OUTPUT_STRIP_TRAILING_WHITESPACE)

string(CONCAT FILE_NAME "_hubbard_I" ${SUFFIX})

file(READ "${CMAKE_CURRENT_BINARY_DIR}/hubbard_I.py" file_content)
# Replace the line
string(REPLACE "import _hubbard_I" "from . import _hubbard_I" modified_content "${file_content}")
# Write the modified content back to the file
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/hubbard_I.py" "${modified_content}")

set(PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/${PROJECT_NAME})
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/${FILE_NAME} DESTINATION ${PYTHON_LIB_DEST})
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/hubbard_I.py DESTINATION ${PYTHON_LIB_DEST})
