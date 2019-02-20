#Set the package name
SET(PACKAGE_NAME "xolotl.modified")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Modified Reaction")
#Set the library name
SET(LIBRARY_NAME "xolotlModified")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS trapmutation/*.h resolution/*.h desorption/*.h)

#Grab all of the source files
file(GLOB SRC trapmutation/*.cpp resolution/*.cpp desorption/*.cpp)

#Include directories
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/modifiedreaction/trapmutation)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/modifiedreaction/resolution)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/modifiedreaction/desorption)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants/neclusters)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/advection)

#Add a library
add_library(${LIBRARY_NAME} STATIC ${SRC})
target_link_libraries(${LIBRARY_NAME} xolotlAdvection)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
