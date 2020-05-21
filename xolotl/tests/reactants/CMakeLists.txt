#Set the package name
SET(PACKAGE_NAME "xolotl.tests.reactants")

#Set the description
SET(PACKAGE_DESCRIPTION "Tests for the Xolotl reactants package")

#Get the test files
file(GLOB tests *Tester.cpp)

#Create the utility library
add_library(testUtils STATIC SimpleReactionNetwork.cpp)
target_link_libraries(testUtils xolotlReactants)

#Make executables and link libraries for testers
foreach(test ${tests})
    message(STATUS "Making test ${test}")
    get_filename_component(testName ${test} NAME_WE)
    add_executable(${testName} ${test})
    target_link_libraries(${testName}
        xolotlReactants
        xolotlCL
        xolotlPerf
        testUtils
        Boost::unit_test_framework
    )
    add_test(${testName} ${testName}) 
    #add a label so the tests can be run separately
    set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})   
endforeach(test ${tests})
