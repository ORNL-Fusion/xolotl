if(VTKm_FOUND)
    file(GLOB tests *Tester.cpp)
else(VTKm_FOUND)
    file(GLOB tests DummyPlotTester.cpp DummyDataProviderTester.cpp)
endif(VTKm_FOUND)

add_tests(tests LIBS xolotlViz LABEL "xolotl.tests.viz")
