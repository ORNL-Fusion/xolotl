file(GLOB tests
    advection/*Tester.cpp
    diffusion/*Tester.cpp
    flux/*Tester.cpp
    modified/*Tester.cpp
    network/*Tester.cpp
    temperature/*Tester.cpp
)

foreach(test ${tests})
    get_filename_component(testName ${test} NAME_WE)
    add_executable(${testName} ${test})
    target_link_libraries(${testName}
        xolotlCore
        xolotlOptions
        xolotlPerf
        Boost::unit_test_framework
    )
    add_test(${testName} ${testName}) 
    #add a label so the tests can be run separately
    set_property(TEST ${testName} PROPERTY LABELS "xolotl.tests.core")
endforeach()
