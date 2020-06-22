set(XOLOTL_INTERFACE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(XOLOTL_INTERFACE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XOLOTL_INTERFACE_HEADER_DIR ${XOLOTL_INTERFACE_INCLUDE_DIR}/xolotl/interface)

set(XOLOTL_INTERFACE_HEADERS
    ${XOLOTL_INTERFACE_HEADER_DIR}/Interface.h
)

set(XOLOTL_INTERFACE_SOURCES
    ${XOLOTL_INTERFACE_SOURCE_DIR}/Interface.cpp
)

add_library(xolotlInterface STATIC
    ${XOLOTL_INTERFACE_SOURCES}
    ${XOLOTL_INTERFACE_HEADERS}
)
target_link_libraries(xolotlInterface PUBLIC
    xolotlCore
    xolotlFactory
    xolotlSolver
    xolotlIO
    xolotlOptions
    xolotlPerf
    xolotlViz
)
target_include_directories(xolotlInterface PUBLIC ${XOLOTL_INTERFACE_INCLUDE_DIR})

install(FILES ${HEADERS} DESTINATION include)