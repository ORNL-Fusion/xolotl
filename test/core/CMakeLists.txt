file(GLOB tests
    advection/*Tester.cpp
    diffusion/*Tester.cpp
    flux/*Tester.cpp
    network/*Tester.cpp
    temperature/*Tester.cpp
)

add_tests(tests LIBS xolotlCore xolotlOptions xolotlPerf LABEL "xolotl.tests.core")
