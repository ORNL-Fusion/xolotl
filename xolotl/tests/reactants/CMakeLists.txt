#Set the package name
SET(PACKAGE_NAME "xolotl.tests.reactants")

#Set the description
SET(PACKAGE_DESCRIPTION "Tests for the Xolotl reactants package")

#include directories from the source and boost binaries
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters)
include_directories(${Boost_BINARY_DIRS})

#link the directories
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlCore/reactants)

#Get the test files
file(GLOB tests *.cpp)

#If boost was found, create tests
if(Boost_FOUND)
    #Make executables and link libraries for testers
    foreach(test ${tests})
        message(STATUS "Making test ${test}")
        get_filename_component(testName ${test} NAME_WE)
        add_executable(${testName} ${test})
        add_test(${testName} ${testName}) 
        #add a label so the tests can be run separately
        set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})   
    endforeach(test ${tests})
endif(Boost_FOUND)