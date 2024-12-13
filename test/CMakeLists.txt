# List of all tests
set(all_tests 
    interact
    solver
    wien2kkpoints
    vaspkpoints
    wien2kshells
    vaspshells
)

foreach(test ${all_tests})
  set (test_dir ${CMAKE_CURRENT_BINARY_DIR}/${test})
  
  file(COPY ${test}/ DESTINATION ${test_dir})
  #foreach(file test.py ref.h5)
  #  FILE(COPY ${test}/${file} DESTINATION ${test_dir})
  #endforeach()

  add_test(NAME ${test}
           #COMMAND bash ${test}.sh
           COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${test_dir}/test.py
           WORKING_DIRECTORY ${test_dir}
          )

  set_property(TEST ${test} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH})
endforeach()
