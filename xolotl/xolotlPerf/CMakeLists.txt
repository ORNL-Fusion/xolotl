#Set the package name
SET(PACKAGE_NAME "xolotl.performance")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Perf")
#Set the library name
SET(LIBRARY_NAME "xolotlPerf")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Collect all of the cpp files in this folder 
#and glob them in SRC
file(GLOB SRC *.cpp)

#Include headers in order to build
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/dummy)

# Find GPTL
FIND_PACKAGE(GPTL)

# Find PAPI
FIND_PACKAGE(PAPI)

if(GPTL_FOUND AND PAPI_FOUND)
    message(STATUS "PAPI includes = ${PAPI_INCLUDE_DIRS}")
    message(STATUS "PAPI libraries = ${PAPI_LIBRARIES}")
    message(STATUS "GPTL includes = ${GPTL_INCLUDE_DIRS}")
    message(STATUS "GPTL libraries = ${GPTL_LIBRARIES}")

    #Add the standard subdirectory
    add_subdirectory(standard)
endif(GPTL_FOUND AND PAPI_FOUND)

#-----------------------------------------------
#if(GPTL_FOUND AND PAPI_FOUND)
#    INCLUDE_DIRECTORIES(${GPTL_INCLUDE_DIRS})
#    INCLUDE_DIRECTORIES(${PAPI_INCLUDE_DIRS})
#    set( hasGPTL "")
#    set( hasPAPI "")
#    message(STATUS "PAPI includes = ${PAPI_INCLUDE_DIRS}")
#    message(STATUS "PAPI libraries = ${PAPI_LIBRARIES}")
#    message(STATUS "GPTL includes = ${GPTL_INCLUDE_DIRS}")
#    message(STATUS "GPTL libraries = ${GPTL_LIBRARIES}")
#    
    #Collect all header filenames in this project 
    #and glob them in HEADERS
#    file(GLOB HEADERS standard/*.h)
#    
    #Collect all of the cpp files in this folder 
    #and glob them in SRC
#    file(GLOB SRC standard/*.cpp)
#endif(GPTL_FOUND AND PAPI_FOUND)
#-------------------------------------------------------

#if(GPTL_FOUND AND PAPI_FOUND)
#    include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/standard)
#    include_directories(${GPTL_INCLUDE_DIRS})
#    include_directories(${PAPI_INCLUDE_DIRS})
#endif(GPTL_FOUND AND PAPI_FOUND)

#Add the library
add_library(${LIBRARY_NAME} SHARED ${SRC})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)
