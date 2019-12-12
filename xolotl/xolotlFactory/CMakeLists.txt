#Set the package name
SET(PACKAGE_NAME "xolotl.factory")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Factory")
#Set the library name
SET(LIBRARY_NAME "xolotlFactory")

#Collect all header filenames in this project
#and glob them in HEADERS
file(GLOB HEADERS material/*.h vizHandler/*.h temperatureHandler/*.h
solverHandler/*.h reactionHandler/*.h)

#Grab all of the source files
file(GLOB SRC material/*.cpp vizHandler/*.cpp temperatureHandler/*.cpp
solverHandler/*.cpp reactionHandler/*.cpp)

#Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_CURRENT_SOURCE_DIR}/material
                    ${CMAKE_SOURCE_DIR}/xolotlCore
                    ${CMAKE_SOURCE_DIR}/xolotlCore/io
                    ${CMAKE_SOURCE_DIR}/xolotlCore/advection
                    ${CMAKE_SOURCE_DIR}/xolotlCore/flux
                    ${CMAKE_SOURCE_DIR}/xolotlCore/diffusion
                    ${CMAKE_SOURCE_DIR}/xolotlCore/commandline
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/neclusters
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/feclusters
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/alloyclusters
                    ${CMAKE_SOURCE_DIR}/xolotlCore/modifiedreaction/trapmutation
                    ${CMAKE_SOURCE_DIR}/xolotlCore/modifiedreaction/resolution
                    ${CMAKE_SOURCE_DIR}/xolotlCore/modifiedreaction/desorption
                    ${CMAKE_SOURCE_DIR}/xolotlCore/modifiedreaction/heterogeneousnucleation
                    ${CMAKE_SOURCE_DIR}/xolotlPerf
                    ${CMAKE_CURRENT_SOURCE_DIR}/vizHandler
                    ${CMAKE_BINARY_DIR}/xolotlFactory
                    ${CMAKE_SOURCE_DIR}/xolotlViz
                    ${CMAKE_SOURCE_DIR}/xolotlViz/dummy
                    ${CMAKE_CURRENT_SOURCE_DIR}/temperatureHandler
                    ${CMAKE_SOURCE_DIR}/xolotlCore/temperature
                    ${CMAKE_CURRENT_SOURCE_DIR}/solverHandler
                    ${CMAKE_SOURCE_DIR}/xolotlSolver
                    ${CMAKE_SOURCE_DIR}/xolotlSolver/solverhandler
                    ${CMAKE_SOURCE_DIR}
                    ${CMAKE_BINARY_DIR}
                    ${Boost_INCLUDE_DIR}
                    ${PETSC_INCLUDES})


# Find VTKm
if(VTKM_DIR)
    include(${CMAKE_SOURCE_DIR}/cmake-modules/SetupVTKm.cmake)
else()
    message(STATUS "Visualization support needs explicit VTKm_DIR.")
endif()

if(VTKM_FOUND)
    include_directories(${CMAKE_SOURCE_DIR}/xolotlViz/standard)
endif(VTKM_FOUND)


#Add a library to hold the factory code
add_library(${LIBRARY_NAME} STATIC ${SRC})
target_link_libraries(${LIBRARY_NAME} xolotlViz xolotlSolver)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
