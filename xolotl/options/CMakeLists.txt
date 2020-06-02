set(XOLOTL_OPTIONS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(XOLOTL_OPTIONS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XOLOTL_OPTIONS_HEADER_DIR ${XOLOTL_OPTIONS_INCLUDE_DIR}/xolotl/options)

set(XOLOTL_OPTIONS_HEADERS
    ${XOLOTL_OPTIONS_HEADER_DIR}/IOptions.h
    ${XOLOTL_OPTIONS_HEADER_DIR}/Options.h
)

set(XOLOTL_OPTIONS_SOURCES
    ${XOLOTL_OPTIONS_SOURCE_DIR}/Options.cpp
)

add_library(xolotlOptions STATIC
    ${XOLOTL_OPTIONS_SOURCES}
    ${XOLOTL_OPTIONS_HEADERS}
)
target_link_libraries(xolotlOptions PUBLIC
    xolotlPerf
    Boost::program_options
    MPI::MPI_CXX
)
target_include_directories(xolotlOptions PUBLIC ${XOLOTL_OPTIONS_INCLUDE_DIR})

install(FILES ${HEADERS} DESTINATION include)
