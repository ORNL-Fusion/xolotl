set(XOLOTL_FACTORY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(XOLOTL_FACTORY_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XOLOTL_FACTORY_INCLUDE_DIR ${XOLOTL_FACTORY_INCLUDE_DIR} PARENT_SCOPE)
set(XOLOTL_FACTORY_HEADER_DIR ${XOLOTL_FACTORY_INCLUDE_DIR}/xolotl/factory)

include(src/material/Include.cmake)
include(src/network/Include.cmake)
include(src/solver/Include.cmake)
include(src/temperature/Include.cmake)
include(src/viz/Include.cmake)

add_library(xolotlFactory STATIC
    ${XOLOTL_FACTORY_SOURCES}
    ${XOLOTL_FACTORY_HEADERS}
)
target_link_libraries(xolotlFactory PUBLIC
    xolotlCore
    xolotlOptions
    xolotlPerf
    xolotlViz
    xolotlSolver
)
target_include_directories(xolotlFactory PUBLIC ${XOLOTL_FACTORY_INCLUDE_DIR})

install(FILES
    ${XOLOTL_FACTORY_HEADER_DIR}/material/IMaterialFactory.h
    DESTINATION include/xolotl/factory/material
)
