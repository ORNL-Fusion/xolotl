#Set the package name
SET(PACKAGE_NAME "xolotl.tests.reactants")

#Set the description
SET(PACKAGE_DESCRIPTION "Tests for the Xolotl reactants package")

#Include directories from the source and boost binaries
include_directories(${CMAKE_SOURCE_DIR}
                    ${CMAKE_SOURCE_DIR}/xolotlCore/io
                    ${CMAKE_SOURCE_DIR}/xolotlCore/commandline
                    ${CMAKE_SOURCE_DIR}/xolotlCore/commandline/optionhandlers
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/neclusters
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/feclusters
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/alloyclusters
                    ${CMAKE_SOURCE_DIR}/xolotlCore
                    ${CMAKE_SOURCE_DIR}/xolotlPerf
                    ${CMAKE_SOURCE_DIR}/xolotlPerf/dummy
                    ${Boost_BINARY_DIRS}
                    ${CMAKE_BINARY_DIR})

#Link the directories
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlCore/reactants)


#Get the test files
file(GLOB tests *Tester.cpp)

#Create the utility library
add_library(testUtils STATIC SimpleReactionNetwork.cpp)
target_link_libraries(testUtils xolotlReactants)

#If boost was found, create tests
if(Boost_FOUND)
    #Make executables and link libraries for testers
    foreach(test ${tests})
        message(STATUS "Making test ${test}")
        get_filename_component(testName ${test} NAME_WE)
        add_executable(${testName} ${test})
        target_link_libraries(${testName} xolotlReactants xolotlCL xolotlPerf testUtils ${Boost_LIBRARIES})
        add_test(${testName} ${testName}) 
        #add a label so the tests can be run separately
        set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})   
    endforeach(test ${tests})
endif(Boost_FOUND)
