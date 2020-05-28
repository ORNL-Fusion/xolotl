if(NOT ${BUILD_TESTING})
    return()
endif()

set(XOLOTL_TEST_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(core)
add_subdirectory(io)
add_subdirectory(options)
add_subdirectory(perf)
add_subdirectory(solver)
add_subdirectory(viz)
add_subdirectory(system)
