#Set the package name
SET(PACKAGE_NAME "xolotl.temperature")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Temperature")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Include directories
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/commandline)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION ${PROJECT_BINARY_DIR}/include)
