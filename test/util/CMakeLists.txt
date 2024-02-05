set(tests
    StepSequenceTester.cpp
    TimeStepperTester.cpp
    TokenizerTester.cpp
)

add_tests(tests LIBS xolotlUtil LABEL "xolotl.tests.util")
