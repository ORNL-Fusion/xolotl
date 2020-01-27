#Set the package name
SET(PACKAGE_NAME "xolotl.tests.commandline")

#Set the description
SET(PACKAGE_DESCRIPTION "Tests for the Xolotl commandline package")

    
# Where to look for include files.
include_directories(${CMAKE_SOURCE_DIR}
                    ${CMAKE_SOURCE_DIR}/xolotlCore
                    ${CMAKE_SOURCE_DIR}/xolotlCore/commandline
                    ${CMAKE_SOURCE_DIR}/xolotlPerf
                    ${Boost_BINARY_DIRS}
                    ${CMAKE_BINARY_DIR}) 

# Figure out which tests to build
file(GLOB tests *Tester.cpp)
    
#If boost was found, create tests
if(Boost_FOUND)

    #Make executables and link libraries for testers
    foreach(test ${tests})
        message(STATUS "Making test ${test}")
        get_filename_component(testName ${test} NAME_WE)
        add_executable(${testName} ${test})
        target_link_libraries(${testName} xolotlCL MPI::MPI_CXX)
        add_test(${testName} ${testName}) 
        #add a label so the tests can be run separately
        set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})   
    endforeach(test ${tests})

endif(Boost_FOUND)

