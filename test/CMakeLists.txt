if(NOT ${BUILD_TESTING})
    return()
endif()

find_package(Boost REQUIRED COMPONENTS unit_test_framework)

set(XOLOTL_TEST_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

get_filename_component(XOLOTL_TEST_DATA_DIR
    "${CMAKE_SOURCE_DIR}/benchmarks"
    ABSOLUTE
)
set(XOLOTL_BUILD_DIR ${CMAKE_BINARY_DIR})

execute_process(COMMAND
    ${CMAKE_COMMAND} -E make_directory ${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/test
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
    "${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/test/config.h"
)
include_directories(${XOLOTL_BINARY_INCLUDE_DIR} ${XOLOTL_TEST_INCLUDE_DIR})

add_subdirectory(util)
add_subdirectory(core)
add_subdirectory(io)
add_subdirectory(interface)
add_subdirectory(options)
add_subdirectory(perf)
add_subdirectory(viz)
add_subdirectory(system)
