set(XOLOTL_CORE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(XOLOTL_CORE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XOLOTL_CORE_HEADER_DIR ${XOLOTL_CORE_INCLUDE_DIR}/xolotl/core)

set(XOLOTL_CORE_HEADERS
    ${XOLOTL_CORE_HEADER_DIR}/Constants.h
    ${XOLOTL_CORE_HEADER_DIR}/InvalidGridDimension.h
)

include(src/advection/Include.cmake)
include(src/diffusion/Include.cmake)
include(src/flux/Include.cmake)
include(src/material/Include.cmake)
include(src/network/Include.cmake)
include(src/temperature/Include.cmake)
include(src/modified/Include.cmake)

add_library(xolotlCore SHARED
    ${XOLOTL_CORE_SOURCES}
    ${XOLOTL_CORE_HEADERS}
    ${XOLOTL_CORE_TPP_FILES}
)
target_link_libraries(xolotlCore PUBLIC
    xolotlFactory
    xolotlPerf
    xolotlOptions
    xolotlUtil
    plsm
    MPI::MPI_CXX
)
target_include_directories(xolotlCore PUBLIC
    $<BUILD_INTERFACE:${XOLOTL_CORE_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>
)
set_target_properties(xolotlCore PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

install(TARGETS xolotlCore EXPORT Xolotl LIBRARY DESTINATION lib)
