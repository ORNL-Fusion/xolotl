add_library(SystemTestCase SystemTestCase.cpp)
target_link_libraries(SystemTestCase PUBLIC
    xolotlIO
    xolotlInterface
    Boost::unit_test_framework
)
target_include_directories(SystemTestCase PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_BINARY_DIR}
)

add_executable(SystemTester SystemTester.cpp)
target_link_libraries(SystemTester SystemTestCase)
add_test(SystemTester SystemTester)
set_property(TEST SystemTester PROPERTY LABELS "xolotl.tests.system")

add_executable(BenchmarkTester BenchmarkTester.cpp)
target_link_libraries(BenchmarkTester SystemTestCase)
