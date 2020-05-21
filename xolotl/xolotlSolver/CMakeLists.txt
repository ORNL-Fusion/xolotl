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

#Add the library
add_library(xolotlSolver STATIC ${SRC})
target_link_libraries(xolotlSolver PUBLIC
    xolotlReactants
    xolotlIO
    xolotlCL
    xolotlDiffusion
    xolotlAdvection
    xolotlFlux
    xolotlModified
    xolotlTemperature
    xolotlFactory
    xolotlPerf
    xolotlViz
    ${PETSC_LIBRARIES}
)
target_include_directories(xolotlSolver PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/solverhandler
    ${CMAKE_CURRENT_SOURCE_DIR}/monitor
    ${PETSC_INCLUDES}
)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS xolotlSolver DESTINATION lib)
