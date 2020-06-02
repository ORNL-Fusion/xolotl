set(XOLOTL_BINARY_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include)
set(XOLOTL_BINARY_INCLUDE_DIR ${XOLOTL_BINARY_INCLUDE_DIR} PARENT_SCOPE)

add_subdirectory(util)
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
