# All Python files. Copy them in the build dir to have a complete package for the tests.
file(GLOB_RECURSE python_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
foreach(file ${python_sources})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# Install python files to proper location
set(PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/${PROJECT_NAME})
#execute_process(COMMAND cp ${CMAKE_CURRENT_BINARY_DIR}/../../fortran/hubbard_I.py ${CMAKE_CURRENT_BINARY_DIR})

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT} FILES_MATCHING PATTERN "*.py" PATTERN "*_desc.py" EXCLUDE)
execute_process(COMMAND ln -fs ${CMAKE_CURRENT_BINARY_DIR}/../../fortran/${FILE_NAME} ${CMAKE_BINARY_DIR}/python/hubbardI)

