set(PACKAGE_NAME "xolotl.tests.system")

set(PACKAGE_DESCRIPTION "Xolotl System Tests")

if(NOT Boost_FOUND)
    return()
endif()

add_library(SystemTest INTERFACE)
target_link_libraries(SystemTest INTERFACE ${Boost_LIBRARIES})
target_include_directories(SystemTest INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/xolotlCore/io
)
get_filename_component(__test_data_dir
    "${CMAKE_SOURCE_DIR}/../benchmarks"
    ABSOLUTE
)
target_compile_definitions(SystemTest INTERFACE
    -DTEST_DATA_DIR="${__test_data_dir}"
    -DBUILD_DIR="${CMAKE_BINARY_DIR}"
)

add_executable(SystemTester SystemTester.cpp)
target_link_libraries(SystemTester SystemTest xolotlIO)
add_test(SystemTester SystemTester)
set_property(TEST SystemTester PROPERTY LABELS ${PACKAGE_NAME})

add_executable(BenchmarkTester BenchmarkTester.cpp)
target_link_libraries(BenchmarkTester SystemTest)
