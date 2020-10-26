set(XOLOTL_SOLVER_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(XOLOTL_SOLVER_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XOLOTL_SOLVER_HEADER_DIR ${XOLOTL_SOLVER_INCLUDE_DIR}/xolotl/solver)

set(XOLOTL_SOLVER_HEADERS
    ${XOLOTL_SOLVER_HEADER_DIR}/ISolver.h
    ${XOLOTL_SOLVER_HEADER_DIR}/Solver.h
    ${XOLOTL_SOLVER_HEADER_DIR}/PetscSolver.h
    ${XOLOTL_SOLVER_HEADER_DIR}/handler/ISolverHandler.h
    ${XOLOTL_SOLVER_HEADER_DIR}/handler/PetscSolver0DHandler.h
    ${XOLOTL_SOLVER_HEADER_DIR}/handler/PetscSolver1DHandler.h
    ${XOLOTL_SOLVER_HEADER_DIR}/handler/PetscSolver2DHandler.h
    ${XOLOTL_SOLVER_HEADER_DIR}/handler/PetscSolver3DHandler.h
    ${XOLOTL_SOLVER_HEADER_DIR}/handler/PetscSolverHandler.h
    ${XOLOTL_SOLVER_HEADER_DIR}/handler/SolverHandler.h
    ${XOLOTL_SOLVER_HEADER_DIR}/monitor/Monitor.h
)

set(XOLOTL_SOLVER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PetscSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Solver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/handler/PetscSolver0DHandler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/handler/PetscSolver1DHandler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/handler/PetscSolver2DHandler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/handler/PetscSolver3DHandler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/handler/PetscSolverHandler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/monitor/Monitor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/monitor/Monitor0D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/monitor/Monitor1D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/monitor/Monitor2D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/monitor/Monitor3D.cpp
)

add_library(xolotlSolver SHARED
    ${XOLOTL_SOLVER_SOURCES}
    ${XOLOTL_SOLVER_HEADERS}
)
target_link_libraries(xolotlSolver PUBLIC
    xolotlCore
    xolotlIO
    xolotlOptions
    xolotlFactory
    xolotlPerf
    xolotlViz
    PETSc
)
target_include_directories(xolotlSolver PUBLIC
    $<BUILD_INTERFACE:${XOLOTL_SOLVER_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>
)
set_target_properties(xolotlSolver PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

install(TARGETS xolotlSolver EXPORT Xolotl LIBRARY DESTINATION lib)
