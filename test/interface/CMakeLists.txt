file(GLOB tests *Tester.cpp)

add_tests(tests LIBS xolotlInterface LABEL "xolotl.tests.interface")
