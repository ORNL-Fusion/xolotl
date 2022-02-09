set(XOLOTL_INTERFACE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(XOLOTL_INTERFACE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XOLOTL_INTERFACE_HEADER_DIR ${XOLOTL_INTERFACE_INCLUDE_DIR}/xolotl/interface)

set(XOLOTL_INTERFACE_HEADERS
    ${XOLOTL_INTERFACE_HEADER_DIR}/Interface.h
)

set(XOLOTL_INTERFACE_SOURCES
    ${XOLOTL_INTERFACE_SOURCE_DIR}/Interface.cpp
)

add_library(xolotlInterface SHARED
    ${XOLOTL_INTERFACE_SOURCES}
    ${XOLOTL_INTERFACE_HEADERS}
)
add_dependencies(xolotlInterface xolotlVersion)
target_link_libraries(xolotlInterface
    PRIVATE
    xolotlCore
    xolotlFactory
    xolotlSolver
    xolotlIO
    xolotlOptions
    xolotlPerf
    xolotlUtil
    xolotlViz
    PUBLIC
    MPI::MPI_CXX
    Kokkos::kokkos
    PETSc
)
target_include_directories(xolotlInterface PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${XOLOTL_INTERFACE_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${XOLOTL_BINARY_INCLUDE_DIR}>
)
set_target_properties(xolotlInterface PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

install(TARGETS xolotlInterface
    EXPORT Xolotl
    LIBRARY DESTINATION lib
)
install(DIRECTORY ${XOLOTL_INTERFACE_INCLUDE_DIR}/xolotl
    DESTINATION include
)
