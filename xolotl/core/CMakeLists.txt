set(XOLOTL_CORE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(XOLOTL_CORE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XOLOTL_CORE_HEADER_DIR ${XOLOTL_CORE_INCLUDE_DIR}/xolotl/core)

set(XOLOTL_CORE_HEADERS
    ${XOLOTL_CORE_HEADER_DIR}/Constants.h
)

include(src/network/Include.cmake)
include(src/diffusion/Include.cmake)
include(src/advection/Include.cmake)
include(src/flux/Include.cmake)
include(src/temperature/Include.cmake)
include(src/modified/Include.cmake)

add_library(xolotlCore STATIC ${XOLOTL_CORE_SOURCES} ${XOLOTL_CORE_HEADERS})
target_link_libraries(xolotlCore PUBLIC
    xolotlPerf
    xolotlOptions
    xolotlUtil
    plsm
    MPI::MPI_CXX
)
target_include_directories(xolotlCore PUBLIC ${XOLOTL_CORE_INCLUDE_DIR})
