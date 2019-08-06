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
include_directories(${CMAKE_SOURCE_DIR}
                    ${CMAKE_SOURCE_DIR}/xolotlSolver
                    ${CMAKE_SOURCE_DIR}/xolotlSolver/solverhandler
                    ${CMAKE_SOURCE_DIR}/xolotlSolver/monitor
                    ${CMAKE_SOURCE_DIR}/xolotlCore
                    ${CMAKE_SOURCE_DIR}/xolotlCore/commandline
                    ${CMAKE_SOURCE_DIR}/xolotlCore/commandline/optionhandlers
                    ${CMAKE_SOURCE_DIR}/xolotlCore/io
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/neclusters
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/feclusters
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants/alloyclusters)
                    ${CMAKE_SOURCE_DIR}/xolotlCore/diffusion
                    ${CMAKE_SOURCE_DIR}/xolotlCore/advection
                    ${CMAKE_SOURCE_DIR}/xolotlCore/flux
                    ${CMAKE_SOURCE_DIR}/xolotlCore/temperature
                    ${CMAKE_SOURCE_DIR}/xolotlCore/modifiedreaction/trapmutation
                    ${CMAKE_SOURCE_DIR}/xolotlCore/modifiedreaction/resolution
                    ${CMAKE_SOURCE_DIR}/xolotlViz
                    ${CMAKE_SOURCE_DIR}/xolotlViz/dataprovider
                    ${CMAKE_SOURCE_DIR}/xolotlPerf
                    ${CMAKE_SOURCE_DIR}/xolotlFactory/material
                    ${CMAKE_SOURCE_DIR}/xolotlFactory/vizHandler
                    ${CMAKE_BINARY_DIR}
                    ${Boost_INCLUDE_DIR}
                    ${PETSC_INCLUDES})

#Add the library
add_library(${LIBRARY_NAME} STATIC ${SRC})
target_link_libraries(${LIBRARY_NAME} xolotlReactants xolotlIO xolotlCL xolotlDiffusion
xolotlAdvection xolotlFlux xolotlModified ${PETSC_LIBRARIES} xolotlPerf xolotlViz)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION ${PROJECT_BINARY_DIR}/include)
