#Set the package name
SET(PACKAGE_NAME "xolotl.flux")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Flux")
#Set the library name
SET(LIBRARY_NAME "xolotlFlux")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Grab all of the source files
file(GLOB SRC *.cpp)

#Include directories
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/commandline)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)

#Add the libraries
add_library(${LIBRARY_NAME} STATIC ${SRC})
target_link_libraries(${LIBRARY_NAME} plsm)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
