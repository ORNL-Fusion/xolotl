#Set the package name
set(PACKAGE_NAME "xolotl.performance")
#Set the description
set(PACKAGE_DESCRIPTION "Xolotl Perf")
#Set the library name
set(LIBRARY_NAME "xolotlPerf")

# Always use the sources and headers common to all registry types.
file(GLOB COMMONHEADERS *.h)
file(GLOB COMMONSRC *.cpp)

# Always include the Dummy (stub) classes.
file(GLOB DUMMYHEADERS *.h dummy/*.h)
file(GLOB DUMMYSRC *.cpp dummy/*.cpp)

# Always include the Standard classes, since we will always be 
# building the support for OS-provided timers at a minimum.
set(STD_HEADERS standard/StdHandlerRegistry.h standard/EventCounter.h)
set(STD_SRC standard/StdHandlerRegistry.cpp)

# Include OS timer support.
# We use the Standard C++ Library interface in <chrono>.
set(OS_HEADERS os/OSHandlerRegistry.h os/OSTimer.h)
set(OS_SRC os/OSHandlerRegistry.cpp os/OSTimer.cpp)

# Check whether PAPI is available.
find_package(PAPI)
if(PAPI_FOUND)
    message(STATUS "PAPI includes = ${PAPI_INCLUDE_DIRS}")
    message(STATUS "PAPI libraries = ${PAPI_LIBRARIES}")

    # FIXME

    include_directories(${PAPI_INCLUDE_DIRS})

    # Add the PAPI-based headers and source files.
    file(GLOB PAPI_HEADERS papi/PAPI*.h)
    file(GLOB PAPI_SRC papi/PAPI*.cpp)

    # Make sure we link with the PAPI library.
    set(MAYBE_PAPI ${PAPI_LIBRARIES})

    # Ensure that our config file indicates we are have the 
    # PAPI-based performance classes available.
    set(HAVE_PAPI 1)

endif(PAPI_FOUND)

set(HEADERS ${COMMONHEADERS} ${DUMMYHEADERS} ${STD_HEADERS} ${OS_HEADERS}
${PAPI_HEADERS})
set(SRC ${COMMONSRC} ${DUMMYSRC} ${STD_SRC} ${OS_SRC} ${PAPI_SRC})

# Specify the library to build
add_library(xolotlPerf STATIC ${SRC})
target_link_libraries(xolotlPerf PUBLIC ${MAYBE_PAPI} MPI::MPI_CXX)
target_include_directories(xolotlPerf PUBLIC
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR}
)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS xolotlPerf DESTINATION lib)

# Configure file for the performance library.
# Note: must do this after all the other checks, or else the 
# contents of the generated file won't take the check results into account.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/perfConfig.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/perfConfig.h"
)
