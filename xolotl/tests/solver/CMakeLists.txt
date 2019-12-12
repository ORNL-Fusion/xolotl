#Set the package name
SET(PACKAGE_NAME "xolotl.tests.solver")

#Set the description
SET(PACKAGE_DESCRIPTION "Tests for the Xolotl solver package")

#Include directories from the source and boost binaries
include_directories(${CMAKE_SOURCE_DIR}
                    ${CMAKE_SOURCE_DIR}/xolotlSolver
                    ${CMAKE_SOURCE_DIR}/xolotlSolver/solverhandler
                    ${CMAKE_SOURCE_DIR}/xolotlCore
                    ${CMAKE_SOURCE_DIR}/xolotlCore/io
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/neclusters
                    ${CMAKE_SOURCE_DIR}/xolotlCore/diffusion
                    ${CMAKE_SOURCE_DIR}/xolotlCore/advection
                    ${CMAKE_SOURCE_DIR}/xolotlCore/flux
                    ${CMAKE_SOURCE_DIR}/xolotlCore/temperature
                    ${CMAKE_SOURCE_DIR}/xolotlCore/commandline
                    ${CMAKE_SOURCE_DIR}/xolotlCore/modifiedreaction/trapmutation
                    ${CMAKE_SOURCE_DIR}/xolotlCore/modifiedreaction/resolution
                    ${CMAKE_SOURCE_DIR}/xolotlCore/modifiedreaction/desorption
                    ${CMAKE_SOURCE_DIR}/xolotlCore/modifiedreaction/heterogeneousnucleation
                    ${CMAKE_SOURCE_DIR}/xolotlPerf
                    ${CMAKE_SOURCE_DIR}/xolotlPerf/dummy
                    ${CMAKE_SOURCE_DIR}/xolotlFactory/material
                    ${CMAKE_SOURCE_DIR}/xolotlFactory/temperatureHandler
                    ${CMAKE_SOURCE_DIR}/xolotlFactory/reactionHandler
                    ${CMAKE_SOURCE_DIR}/xolotlFactory/vizHandler
                    ${CMAKE_SOURCE_DIR}/xolotlViz
                    ${Boost_BINARY_DIRS}
                    ${PETSC_INCLUDES}
                    ${CMAKE_BINARY_DIR})

#Link the directories
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlSolver)

#Get the test files
file(GLOB tests *Tester.cpp)

#If boost was found, create tests
if(Boost_FOUND)
    #Make executables and link libraries for testers
    foreach(test ${tests})
        message(STATUS "Making test ${test}")
        get_filename_component(testName ${test} NAME_WE)
        add_executable(${testName} ${test})
        target_link_libraries(${testName} xolotlSolver xolotlFactory 
        ${PETSC_LIBRARIES} ${Boost_LIBRARIES})
        add_test(${testName} ${testName}) 
        #add a label so the tests can be run separately
        set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})   
    endforeach(test ${tests})
endif(Boost_FOUND)

