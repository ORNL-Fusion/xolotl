#Set the package name
SET(PACKAGE_NAME "xolotl.viz")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Visualization")
#Set the library name
SET(LIBRARY_NAME "xolotlViz")

#Include headers in order to build
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/dataprovider)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/labelprovider)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/dummy)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Find EAVL
FIND_PACKAGE(EAVL)

# Find MESA
FIND_PACKAGE(MESA)

# Always include the Dummy (stub) classes.
file(GLOB HEADERS *.h dummy/*.h dataprovider/*.h labelprovider/*.h)
file(GLOB SRC *.cpp dummy/*.cpp dataprovider/*.cpp)

if(EAVL_FOUND AND MESA_FOUND)
    message(STATUS "EAVL includes = ${EAVL_INCLUDE_DIRS}")
    message(STATUS "EAVL config = ${EAVL_CONFIG_DIRS}")
    message(STATUS "EAVL libraries = ${EAVL_LIBRARIES}")
    message(STATUS "Mesa includes = ${MESA_INCLUDE_DIRS}")
    message(STATUS "Mesa libraries = ${MESA_LIBRARIES}")
    
    
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/standard)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/standard/plot)
    include_directories(${EAVL_INCLUDE_DIRS})
    include_directories(${EAVL_INCLUDE_DIRS}/common)
    include_directories(${EAVL_INCLUDE_DIRS}/math)
    include_directories(${EAVL_INCLUDE_DIRS}/rendering)
    include_directories(${EAVL_INCLUDE_DIRS}/operations)
    include_directories(${EAVL_INCLUDE_DIRS}/fonts)
    include_directories(${EAVL_INCLUDE_DIRS}/importers)
    include_directories(${EAVL_CONFIG_DIRS})
    include_directories(${MESA_INCLUDE_DIRS})
    
    # Add the standard classes to the HEADERS and SRCS.
    file(GLOB STDHEADERS standard/*.h standard/plot/*.h)
    set(HEADERS ${HEADERS} ${STDHEADERS})
    file(GLOB STDSRC standard/*.cpp standard/plot/*.cpp)
    set(SRC ${SRC} ${STDSRC})
    
    set(MAYBE_EAVL_MESA ${EAVL_LIBRARIES} ${MESA_LIBRARIES})

    # Ensure that our config file indicates we have built the 
    # standard visualization classes in our visualization infrastructure.
    set(HAVE_VIZLIB_STD 1)
    
endif(EAVL_FOUND AND MESA_FOUND)

#Add the library
add_library(${LIBRARY_NAME} STATIC ${SRC})
message('${SRC}')
message('${HEADERS}')
target_link_libraries(${LIBRARY_NAME} ${MAYBE_EAVL_MESA})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)

# Configure file for the performance library.
# Note: must do this after all the other checks, or else the 
# contents of the generated file won't take the check results into account.
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/XolotlConfigViz.h.in" "${CMAKE_CURRENT_BINARY_DIR}/XolotlConfigViz.h" )



