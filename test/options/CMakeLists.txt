list(APPEND tests
    OptionsTester.cpp
    ConfOptionsTester.cpp
    JSONOptionsTester.cpp
)

add_tests(tests LIBS xolotlOptions LABEL "xolotl.tests.options")
