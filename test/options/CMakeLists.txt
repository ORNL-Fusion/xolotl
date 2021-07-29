file(GLOB tests *Tester.cpp)
    
add_tests(tests LIBS xolotlOptions LABEL "xolotl.tests.options")
