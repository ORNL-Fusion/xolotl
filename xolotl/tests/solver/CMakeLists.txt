#Set the package name
SET(PACKAGE_NAME "xolotl.tests.solver")

#Set the description
SET(PACKAGE_DESCRIPTION "Tests for the Xolotl solver package")

#Include directories from the source and boost binaries
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/xolotlSolver)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/diffusion)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/advection)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/flux)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/temperature)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/commandline)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/dummy)
include_directories(${CMAKE_SOURCE_DIR}/xolotlFactory/material)
include_directories(${Boost_BINARY_DIRS})
include_directories(${PETSC_INCLUDES})

#Link the directories
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlSolver)
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlSolver/flux)
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlSolver/temperature)

#Get the test files
file(GLOB tests *Tester.cpp)

#If boost was found, create tests
if(Boost_FOUND)
    #Make executables and link libraries for testers
    foreach(test ${tests})
        message(STATUS "Making test ${test}")
        get_filename_component(testName ${test} NAME_WE)
        add_executable(${testName} ${test})
        target_link_libraries(${testName} xolotlSolver xolotlReactants ${PETSC_LIBRARIES} xolotlPerf)
        add_test(${testName} ${testName}) 
        #add a label so the tests can be run separately
        set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})   
    endforeach(test ${tests})
endif(Boost_FOUND)

