add_library(SystemTest INTERFACE)
target_link_libraries(SystemTest INTERFACE
    xolotlIO
    Boost::unit_test_framework
)
target_include_directories(SystemTest INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/xolotlCore/io
)
add_dependencies(SystemTest xolotl)

add_executable(SystemTester SystemTester.cpp)
target_link_libraries(SystemTester SystemTest xolotlIO)
add_test(SystemTester SystemTester)
set_property(TEST SystemTester PROPERTY LABELS "xolotl.tests.system")

add_executable(BenchmarkTester BenchmarkTester.cpp)
target_link_libraries(BenchmarkTester SystemTest)
