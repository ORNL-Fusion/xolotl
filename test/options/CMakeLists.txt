list(APPEND tests
    OptionsTester.cpp
)

if(Xolotl_ENABLE_JSON)
    list(APPEND tests
        JSONOptionsTester.cpp
    )
endif()
    
add_tests(tests LIBS xolotlOptions LABEL "xolotl.tests.options")
