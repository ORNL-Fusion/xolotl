#Set the package name
SET(PACKAGE_NAME "xolotl.solver")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Solver")
#Set the library name
SET(LIBRARY_NAME "xolotlSolver")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Collect all of the cpp files in this folder 
#and glob them in SRC
file(GLOB SRC *.cpp)

#Include headers so that the solvers can be built
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters)
include_directories(${PETSC_INCLUDES})

#Add the library
add_library(${LIBRARY_NAME} SHARED ${SRC})
target_link_libraries(${LIBRARY_NAME} xolotlReactants xolotlIO ${PETSC_LIBRARIES} xolotlPerf)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)
