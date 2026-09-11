# All Python files. Copy them in the build dir to have a complete package for the tests.
file(GLOB_RECURSE python_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
foreach(file ${python_sources})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# Install python files to proper location
set(PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/${PROJECT_NAME})
#execute_process(COMMAND cp ${CMAKE_CURRENT_BINARY_DIR}/../../fortran/hubbard_I.py ${CMAKE_CURRENT_BINARY_DIR})

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT} FILES_MATCHING PATTERN "*.py" PATTERN "*_desc.py" EXCLUDE)
# The python sources import the Fortran extension as `hubbardI.hubbard_I`, so
# expose the fortran build directory under that name inside the build tree.
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                  ${CMAKE_BINARY_DIR}/fortran/MagInt
                  ${CMAKE_BINARY_DIR}/python/hubbardI
                COMMAND_ERROR_IS_FATAL ANY)

