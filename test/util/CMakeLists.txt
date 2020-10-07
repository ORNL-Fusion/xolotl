file(GLOB tests *Tester.cpp)

foreach(test ${tests})
    get_filename_component(testName ${test} NAME_WE)
    add_executable(${testName} ${test})
    target_link_libraries(${testName}
        xolotlUtil
        Boost::unit_test_framework
    )
    add_test(${testName} ${testName}) 
    #add a label so the tests can be run separately
    set_property(TEST ${testName} PROPERTY LABELS "xolotl.tests.util")
endforeach()
