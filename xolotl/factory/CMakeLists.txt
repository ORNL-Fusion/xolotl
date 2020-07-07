set(XOLOTL_FACTORY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(XOLOTL_FACTORY_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XOLOTL_FACTORY_HEADER_DIR ${XOLOTL_FACTORY_INCLUDE_DIR}/xolotl/factory)

include(src/material/Include.cmake)
include(src/network/Include.cmake)
include(src/solver/Include.cmake)
include(src/temperature/Include.cmake)
include(src/viz/Include.cmake)

add_library(xolotlFactory SHARED
    ${XOLOTL_FACTORY_SOURCES}
    ${XOLOTL_FACTORY_HEADERS}
)
target_link_libraries(xolotlFactory PUBLIC
    xolotlOptions
)
target_include_directories(xolotlFactory PUBLIC
    $<BUILD_INTERFACE:${XOLOTL_FACTORY_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>
)
set_target_properties(xolotlFactory PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
)

install(TARGETS xolotlFactory EXPORT Xolotl LIBRARY DESTINATION lib)
