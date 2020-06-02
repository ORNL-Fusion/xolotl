file(GLOB tests
    *Tester.cpp
    dummy/*Tester.cpp
    os/*Tester.cpp
    standard/*Tester.cpp
)

if(PAPI_FOUND)
    file(GLOB PAPI_TEST_SRCS papi/*Tester.cpp)
    list(APPEND tests ${PAPI_TEST_SRCS})
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

foreach(test ${tests})
    get_filename_component(testName ${test} NAME_WE)
    add_executable(${testName} ${test})
    target_link_libraries(${testName} xolotlPerf Boost::unit_test_framework)
    add_test(${testName} ${testName})
    #add a label so the tests can be run separately
    set_property(TEST ${testName} PROPERTY LABELS "xolotl.tests.perf")
endforeach(test ${tests})
