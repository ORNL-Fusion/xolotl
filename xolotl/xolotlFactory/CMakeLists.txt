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

# Find VTKm
if(VTKM_DIR)
    include(${CMAKE_SOURCE_DIR}/cmake-modules/SetupVTKm.cmake)
else()
    message(STATUS "Visualization support needs explicit VTKm_DIR.")
endif()

#Add a library to hold the factory code
add_library(xolotlFactory STATIC ${SRC})
target_link_libraries(xolotlFactory PUBLIC
    xolotlIO
    xolotlAdvection
    xolotlFlux
    xolotlDiffusion
    xolotlCL
    xolotlReactants
    xolotlModified
    xolotlPerf
    xolotlViz
    xolotlSolver
)
target_include_directories(xolotlFactory PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/material
    ${CMAKE_CURRENT_SOURCE_DIR}/reactionHandler
    ${CMAKE_CURRENT_SOURCE_DIR}/solverHandler
    ${CMAKE_CURRENT_SOURCE_DIR}/temperatureHandler
    ${CMAKE_CURRENT_SOURCE_DIR}/vizHandler
    ${CMAKE_CURRENT_BINARY_DIR}
)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
