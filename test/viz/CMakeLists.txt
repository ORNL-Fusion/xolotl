if(VTKM_FOUND)
    file(GLOB tests *Tester.cpp)
else(VTKM_FOUND)
    file(GLOB tests DummyPlotTester.cpp DummyDataProviderTester.cpp)
endif(VTKM_FOUND)

foreach(test ${tests})
    get_filename_component(testName ${test} NAME_WE)
    add_executable(${testName} ${test})
    target_link_libraries(${testName}
        xolotlViz
        Boost::unit_test_framework
    )
    add_test(${testName} ${testName})
    #add a label so the tests can be run separately
    set_property(TEST ${testName} PROPERTY LABELS "xolotl.tests.viz")
endforeach()
