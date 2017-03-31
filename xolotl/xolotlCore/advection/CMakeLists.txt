#Set the package name
SET(PACKAGE_NAME "xolotl.advection")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Advection")
#Set the library name
SET(LIBRARY_NAME "xolotlAdvection")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Grab all of the source files
file(GLOB SRC *.cpp)

#Include directories
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters)

#Add a library
add_library(${LIBRARY_NAME} STATIC ${SRC})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
