if(BUILD_DOC)
    # Find the Sphinx executable
    find_program(SPHINX_EXECUTABLE sphinx-build)

    if(NOT SPHINX_EXECUTABLE)
        message(FATAL_ERROR "Sphinx executable not found. Please install Sphinx or specify its location.")
    endif()

    # Define input and output directories for Sphinx
    set(SPHINX_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)
    set(SPHINX_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/html)

    # Point 4: Gather .rst files as dependencies
    file(GLOB_RECURSE SPHINX_RST_FILES ${SPHINX_SOURCE_DIR}/*.rst)

    # Add custom target for building documentation
    add_custom_target(
        doc ALL  # Makes this target part of the default build
        COMMAND ${SPHINX_EXECUTABLE} -b html ${SPHINX_SOURCE_DIR} ${SPHINX_BUILD_DIR}
        DEPENDS ${SPHINX_RST_FILES}  # Rebuild when any .rst file changes
        COMMENT "Generating documentation with Sphinx"
    )

    # Point 5: Add installation rule for documentation
    install(
        DIRECTORY ${SPHINX_BUILD_DIR}/ DESTINATION share/doc/myproject
        COMPONENT documentation
        FILES_MATCHING PATTERN "*.html"  # Only install HTML files
    )
endif()
