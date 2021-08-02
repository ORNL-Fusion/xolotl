list(APPEND tests
    EventCounterTester.cpp
    PerfHandlerTester.cpp
    dummy/DummyEventCounterTester.cpp
    dummy/DummyHardwareCounterTester.cpp
    dummy/DummyTimerTester.cpp
    os/OSTimerTester.cpp
)

if(PAPI_FOUND)
    list(APPEND tests
        papi/PAPIHardwareCounterTester.cpp
        papi/PAPITimerTester.cpp
    )
endif(PAPI_FOUND)

if(CMAKE_BUILD_TYPE MATCHES "^Debug$")
    set(XOLOTL_TEST_HWCTR_DEBUGEXP 1)
    set(XOLOTL_TEST_HWCTR_EXPVAL 1150)
    set(XOLOTL_TEST_HWCTR_TOL 15)
else(CMAKE_BUILD_TYPE MATCHES "^Debug$")
    set(XOLOTL_TEST_HWCTR_DEBUGEXP 0)
    set(XOLOTL_TEST_HWCTR_EXPVAL 1000)
    set(XOLOTL_TEST_HWCTR_TOL 5)
endif(CMAKE_BUILD_TYPE MATCHES "^Debug$")

execute_process(COMMAND
    ${CMAKE_COMMAND} -E make_directory
    ${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/test/perf
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
    "${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/test/perf/config.h"
)

add_tests(tests LIBS xolotlPerf LABEL "xolotl.tests.perf")
