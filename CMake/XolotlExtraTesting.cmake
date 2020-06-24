# Code Coverage
option(CODE_COVERAGE "Add build flags for code coverage?" OFF)
if(${CODE_COVERAGE})
    set(TEST_COVERAGE_OUTPUT_DIR "${CMAKE_BINARY_DIR}/test/coverage")
    execute_process(COMMAND
        ${CMAKE_COMMAND} -E make_directory ${TEST_COVERAGE_OUTPUT_DIR}
    )
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
            message(WARNING "Code coverage results with an optimized \
                (non-Debug) build may be misleading")
        endif()

        # Add flags so that `ctest -T Coverage` will work
        add_compile_options(--coverage)
        add_link_options(--coverage)

        # Add custom target to produce html results using lcov and genhtml
        find_program(LCOV_EXECUTABLE lcov)
        find_program(GENHTML_EXECUTABLE genhtml)
        if(LCOV_EXECUTABLE AND GENHTML_EXECUTABLE)
            set(LCOV_INFO_FILE "${TEST_COVERAGE_OUTPUT_DIR}/testCoverage.info")
            add_custom_target(testCoverage
                COMMAND
                ${LCOV_EXECUTABLE} --capture --no-external
                    --directory ${CMAKE_SOURCE_DIR}
                    --directory ${CMAKE_BINARY_DIR}
                    --output-file ${LCOV_INFO_FILE}
                COMMAND
                ${GENHTML_EXECUTABLE} -o ${TEST_COVERAGE_OUTPUT_DIR}/html
                    ${LCOV_INFO_FILE}
            )
        endif()
    endif()
endif()


