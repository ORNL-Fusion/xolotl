#Set the package name
SET(PACKAGE_NAME "xolotl.tests.advection")

#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Advection tests")

#include directories from the source and boost binaries
include_directories(${CMAKE_SOURCE_DIR}
                    ${CMAKE_SOURCE_DIR}/xolotlCore
                    ${CMAKE_SOURCE_DIR}/xolotlCore/commandline
                    ${CMAKE_SOURCE_DIR}/xolotlCore/commandline/optionhandlers
                    ${CMAKE_SOURCE_DIR}/xolotlCore/advection
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters
                    ${CMAKE_SOURCE_DIR}/xolotlPerf
                    ${CMAKE_SOURCE_DIR}/xolotlPerf/dummy
                    ${Boost_BINARY_DIRS}
                    ${CMAKE_BINARY_DIR})

#link the directories
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlCore/advection)

#Get the test files
file(GLOB tests *.cpp)

#If boost was found, create tests
if(Boost_FOUND)
    #Make executables and link libraries for testers
    foreach(test ${tests})
        message(STATUS "Making test ${test}")
        get_filename_component(testName ${test} NAME_WE)
        add_executable(${testName} ${test})
        target_link_libraries(${testName} xolotlAdvection xolotlPerf xolotlCL
        xolotlReactants ${Boost_LIBRARIES})
        add_test(${testName} ${testName}) 
        #add a label so the tests can be run separately
        set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})   
    endforeach(test ${tests})
endif(Boost_FOUND)
