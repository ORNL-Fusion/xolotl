#Set the package name
SET(PACKAGE_NAME "xolotl.performance")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Perf")
#Set the library name
SET(LIBRARY_NAME "xolotlPerf")

#Include headers in order to build
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_BINARY_DIR})

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

# Check for OS timer support...
include(CheckIncludeFiles)
include(CheckSymbolExists)

# ...first, check for gettimeofday...
check_include_files(sys/time.h HAVE_SYS_TIME_H)
if(HAVE_SYS_TIME_H)
    check_symbol_exists(gettimeofday sys/time.h HAVE_GETTIMEOFDAY)
    if(HAVE_GETTIMEOFDAY)
        set(OS_TIMESTAMP_HDRS ${OS_TIMESTAMP_HDRS} "os/GTODTimestamp.h")
    endif(HAVE_GETTIMEOFDAY)
endif(HAVE_SYS_TIME_H)

# ...next, check for clock_gettime.
check_include_files(time.h HAVE_TIME_H)
if(HAVE_TIME_H)
    set(SAVED_CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
    set(CMAKE_REQUIRED_LIBRARIES ${SAVED_CMAKE_REQUIRED_LIBRARIES} -lrt)
    check_symbol_exists(clock_gettime time.h HAVE_CLOCK_GETTIME)
    if(HAVE_CLOCK_GETTIME)
        set(MAYBE_RT -lrt)
        set(OS_TIMESTAMP_HDRS ${OS_TIMESTAMP_HDRS} "os/CGTTimstamp.h")
    endif(HAVE_CLOCK_GETTIME)
    set(CMAKE_REQUIRED_LIBRARIES ${SAVED_CMAKE_REQUIRED_LIBRARIES})
endif(HAVE_TIME_H)

set(OS_HEADERS os/OSHandlerRegistry.h os/OSTimer.h ${OS_TIMESTAMP_HDRS})
set(OS_SRC os/OSHandlerRegistry.cpp os/OSTimer.cpp)

# Check whether PAPI is available.
FIND_PACKAGE(PAPI)
if(PAPI_FOUND)
    message(STATUS "PAPI includes = ${PAPI_INCLUDE_DIRS}")
    message(STATUS "PAPI libraries = ${PAPI_LIBRARIES}")

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
add_library(${LIBRARY_NAME} STATIC ${SRC})
target_link_libraries(${LIBRARY_NAME} ${MAYBE_PAPI} ${MAYBE_RT})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)


# Configure file for the performance library.
# Note: must do this after all the other checks, or else the 
# contents of the generated file won't take the check results into account.
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/perfConfig.h.in" "${CMAKE_CURRENT_BINARY_DIR}/perfConfig.h")


