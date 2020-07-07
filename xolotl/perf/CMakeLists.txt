if(PAPI_FOUND)
    # Ensure that our config file indicates we are have the 
    # PAPI-based performance classes available.
    set(HAVE_PAPI 1)
endif(PAPI_FOUND)
execute_process(COMMAND
    ${CMAKE_COMMAND} -E make_directory ${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/perf
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
    "${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/perf/config.h"
)

set(XOLOTL_PERF_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(XOLOTL_PERF_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XOLOTL_PERF_HEADER_DIR ${XOLOTL_PERF_INCLUDE_DIR}/xolotl/perf)

set(XOLOTL_PERF_HEADERS
    ${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/perf/config.h
    ${XOLOTL_PERF_HEADER_DIR}/xolotlPerf.h
    ${XOLOTL_PERF_HEADER_DIR}/EventCounter.h
    ${XOLOTL_PERF_HEADER_DIR}/IEventCounter.h
    ${XOLOTL_PERF_HEADER_DIR}/IHandlerRegistry.h
    ${XOLOTL_PERF_HEADER_DIR}/IHardwareCounter.h
    ${XOLOTL_PERF_HEADER_DIR}/ITimer.h
    ${XOLOTL_PERF_HEADER_DIR}/PerfObjStatistics.h
    ${XOLOTL_PERF_HEADER_DIR}/RuntimeError.h
)

set(XOLOTL_PERF_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xolotlPerf.cpp
)

include(src/dummy/Include.cmake)
include(src/os/Include.cmake)
include(src/standard/Include.cmake)
if(PAPI_FOUND)
    include(src/papi/Include.cmake)
endif(PAPI_FOUND)

add_library(xolotlPerf SHARED
    ${XOLOTL_PERF_SOURCES}
    ${XOLOTL_PERF_HEADERS}
)
target_link_libraries(xolotlPerf PUBLIC
    xolotlUtil
)
target_include_directories(xolotlPerf PUBLIC
    $<BUILD_INTERFACE:${XOLOTL_PERF_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${XOLOTL_BINARY_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>
)
set_target_properties(xolotlPerf PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

if(PAPI_FOUND)
    target_link_libraries(xolotlPerf PUBLIC ${PAPI_LIBRARIES})
    target_include_directories(xolotlPerf PUBLIC ${PAPI_INCLUDE_DIRS})
endif(PAPI_FOUND)

install(TARGETS xolotlPerf EXPORT Xolotl LIBRARY DESTINATION lib)
