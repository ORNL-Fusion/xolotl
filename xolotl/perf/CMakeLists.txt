set(XOLOTL_PERF_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(XOLOTL_PERF_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XOLOTL_PERF_HEADER_DIR ${XOLOTL_PERF_INCLUDE_DIR}/xolotl/perf)

set(XOLOTL_PERF_HEADERS
    ${XOLOTL_PERF_HEADER_DIR}/EventCounter.h
    ${XOLOTL_PERF_HEADER_DIR}/IEventCounter.h
    ${XOLOTL_PERF_HEADER_DIR}/IHardwareCounter.h
    ${XOLOTL_PERF_HEADER_DIR}/IPerfHandler.h
    ${XOLOTL_PERF_HEADER_DIR}/ITimer.h
    ${XOLOTL_PERF_HEADER_DIR}/PerfHandler.h
    ${XOLOTL_PERF_HEADER_DIR}/PerfObjStatistics.h
    ${XOLOTL_PERF_HEADER_DIR}/RuntimeError.h
    ${XOLOTL_PERF_HEADER_DIR}/ScopedTimer.h
)

set(XOLOTL_PERF_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PerfHandler.cpp
)

include(src/dummy/Include.cmake)
include(src/os/Include.cmake)
if(PAPI_FOUND)
    include(src/papi/Include.cmake)
endif(PAPI_FOUND)

add_library(xolotlPerf SHARED
    ${XOLOTL_PERF_SOURCES}
    ${XOLOTL_PERF_HEADERS}
)
target_link_libraries(xolotlPerf PUBLIC
    xolotlFactory
    xolotlUtil
)
target_include_directories(xolotlPerf PUBLIC
    $<BUILD_INTERFACE:${XOLOTL_PERF_INCLUDE_DIR}>
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
