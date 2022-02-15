set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(XOLOTL_BINARY_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include)
set(XOLOTL_BINARY_INCLUDE_DIR ${XOLOTL_BINARY_INCLUDE_DIR} PARENT_SCOPE)

configure_file(${CMAKE_CURRENT_LIST_DIR}/config.h.in
    ${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/config.h
)
configure_file(${CMAKE_CURRENT_LIST_DIR}/version.h.in
    ${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/version.h
)
add_custom_target(xolotlVersion ALL
    DEPENDS "${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/version.h"
    SOURCES "${CMAKE_CURRENT_LIST_DIR}/version.h.in"
)

add_subdirectory(util)
add_subdirectory(perf)
add_subdirectory(core)
add_subdirectory(io)
add_subdirectory(options)
add_subdirectory(viz)
add_subdirectory(factory)
add_subdirectory(solver)
add_subdirectory(interface)

add_executable(xolotl main.cpp)
target_link_libraries(xolotl
    xolotlInterface
    Kokkos::kokkos
)
set_target_properties(xolotl PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
)

install(FILES
    ${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/config.h
    DESTINATION include/xolotl
)

install(TARGETS xolotl
    EXPORT Xolotl
    RUNTIME DESTINATION bin
)
install(EXPORT Xolotl DESTINATION lib/Xolotl)

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_SOURCE_DIR}/CMake/XolotlConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/CMake/XolotlConfig.cmake
    INSTALL_DESTINATION lib/Xolotl
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/CMake/XolotlConfig.cmake
    ${CMAKE_SOURCE_DIR}/CMake/FindPETSc.cmake
    DESTINATION lib/Xolotl
)
