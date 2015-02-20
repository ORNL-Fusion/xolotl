#Set the package name
SET(PACKAGE_NAME "xolotl.solver")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Solver")
#Set the library name
SET(LIBRARY_NAME "xolotlSolver")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h solverhandler/*.h monitor/*.h)

#Collect all of the cpp files in this folder 
#and glob them in SRC
file(GLOB SRC *.cpp solverhandler/*.cpp monitor/*.cpp)

#Include headers so that the solvers can be built
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/xolotlSolver)
include_directories(${CMAKE_SOURCE_DIR}/xolotlSolver/solverhandler)
include_directories(${CMAKE_SOURCE_DIR}/xolotlSolver/monitor)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/commandline)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/io)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/diffusion)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/advection)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/flux)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/temperature)
include_directories(${CMAKE_SOURCE_DIR}/xolotlViz)
include_directories(${CMAKE_SOURCE_DIR}/xolotlViz/dataprovider)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)
include_directories(${CMAKE_SOURCE_DIR}/xolotlFactory/material)
include_directories(${CMAKE_SOURCE_DIR}/xolotlFactory/vizHandler)

include_directories(${PETSC_INCLUDES})

#Add the library
add_library(${LIBRARY_NAME} STATIC ${SRC})
target_link_libraries(${LIBRARY_NAME} xolotlReactants xolotlIO xolotlCL xolotlDiffusion
xolotlAdvection xolotlFlux xolotlTemperature ${PETSC_LIBRARIES} xolotlPerf xolotlViz)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)
