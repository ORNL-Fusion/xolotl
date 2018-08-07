#Set the package name
SET(PACKAGE_NAME "xolotl.reactants")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Reactants")
#Set the library name
SET(LIBRARY_NAME "xolotlReactants")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h psiclusters/*.h neclusters/*.h feclusters/*.h)

#Collect all of the cpp files in this folder 
#and the psiclusters folder and glob them in SRC
file(GLOB SRC *.cpp psiclusters/*.cpp neclusters/*.cpp feclusters/*.cpp)

#Include this directory so that the psiclusters can be built
include_directories(${CMAKE_SOURCE_DIR}
                    ${CMAKE_SOURCE_DIR}/xolotlCore
                    ${CMAKE_SOURCE_DIR}/xolotlCore/io
                    ${CMAKE_SOURCE_DIR}/xolotlCore/commandline
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/neclusters
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/feclusters
                    ${CMAKE_SOURCE_DIR}/xolotlPerf
                    ${CMAKE_SOURCE_DIR}/xolotlPerf/dummy
                    ${CMAKE_SOURCE_DIR}/xolotlPerf/standard
                    ${CMAKE_BINARY_DIR}
                    ${Boost_INCLUDE_DIR})

#Add the library
add_library(${LIBRARY_NAME} STATIC ${SRC})
target_link_libraries(${LIBRARY_NAME} xolotlPerf xolotlIO)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)
