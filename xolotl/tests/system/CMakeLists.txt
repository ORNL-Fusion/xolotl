set(PACKAGE_NAME "xolotl.tests.system")

set(PACKAGE_DESCRIPTION "Xolotl System Tests")

if(NOT Boost_FOUND)
    return()
endif()

add_executable(SystemTester SystemTester.cpp)
target_link_libraries(SystemTester ${Boost_LIBRARIES})
get_filename_component(__test_data_dir
    "${CMAKE_SOURCE_DIR}/../benchmarks"
    ABSOLUTE
)
target_compile_definitions(SystemTester PUBLIC
    -DTEST_DATA_DIR="${__test_data_dir}"
    -DBUILD_DIR="${CMAKE_BINARY_DIR}"
)
