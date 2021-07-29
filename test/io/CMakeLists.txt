file(GLOB tests *Tester.cpp)

add_tests(tests LIBS xolotlIO LABEL "xolotl.tests.io")
