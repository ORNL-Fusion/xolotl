#Set the package name
SET(PACKAGE_NAME "xolotl.reactants")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Reactants")
#Set the library name
SET(LIBRARY_NAME "xolotlReactants")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h psiclusters/*.h)

#Collect all of the cpp files in this folder 
#and the psiclusters folder and glob them in SRC
file(GLOB SRC *.cpp psiclusters/*.cpp)

#Include this directory so that the psiclusters can be built
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/io)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/dummy)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/standard)

#Add the library
add_library(${LIBRARY_NAME} STATIC ${SRC})
target_link_libraries(${LIBRARY_NAME} xolotlPerf)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)
