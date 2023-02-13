list(APPEND tests
    DummyPlotTester.cpp
)
if(VTKm_FOUND)
    list(APPEND tests
        CvsXDPTester.cpp
        CvsXYDPTester.cpp
        DataProviderTester.cpp
        PlotTester.cpp
    )
endif(VTKm_FOUND)

add_tests(tests LIBS xolotlViz LABEL "xolotl.tests.viz")
