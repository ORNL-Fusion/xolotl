set(XOLOTL_BINARY_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include)
execute_process(COMMAND
    ${CMAKE_COMMAND} -E make_directory ${XOLOTL_BINARY_INCLUDE_DIR}/xolotl
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
    "${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/config.h"
)

add_subdirectory(perf)
add_subdirectory(core)
add_subdirectory(io)
add_subdirectory(options)
add_subdirectory(viz)
add_subdirectory(factory)
add_subdirectory(solver)

add_executable(xolotl main.cpp)
target_link_libraries(xolotl
    xolotlCore
    xolotlFactory
    xolotlSolver
    xolotlIO
    xolotlOptions
    xolotlPerf
    xolotlViz
)
target_include_directories(xolotl PUBLIC ${XOLOTL_BINARY_INCLUDE_DIR})
