#Set the package name
SET(PACKAGE_NAME "xolotl.tests.temperature")

#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Temperature tests")

#Get the test files
file(GLOB tests *.cpp)

#Make executables and link libraries for testers
foreach(test ${tests})
    message(STATUS "Making test ${test}")
    get_filename_component(testName ${test} NAME_WE)
    add_executable(${testName} ${test})
    target_link_libraries(${testName}
        xolotlReactants
        xolotlCL
        xolotlTemperature
        Boost::unit_test_framework
    )
    add_test(${testName} ${testName}) 
    #add a label so the tests can be run separately
    set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})   
endforeach(test ${tests})
