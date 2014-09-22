#Set the package name
SET(PACKAGE_NAME "xolotl.tests.performance")

#Set the description
SET(PACKAGE_DESCRIPTION "Tests for the Xolotl performance package")

# Our tests depend on Boost.
if(Boost_FOUND)

    # We will need the Boost headers.
    include_directories(${Boost_BINARY_DIRS})

    # Look for headers in the source tree and binary tree.
    include_directories(${CMAKE_SOURCE_DIR})
    include_directories(${CMAKE_BINARY_DIR})

    # Always build tests for the dummy classes.
    file(GLOB DUMMY_TEST_SRCS Dummy*Tester.cpp)

    # Always build the testers for the Standard classes that are always built
    set(COMMON_TEST_SRCS EventCounterTester.cpp StdHandlerRegistryTester.cpp)

    # Always build the testers for the OS classes that are always built.
    file(GLOB OS_TEST_SRCS OS*Tester.cpp)

    # Tests for the standard classes.
    # We re-do the test for PAPI here because cmake's scoping rules
    # mean that PAPI_FOUND (set earlier in the xolotlPerf library's
    # CMakeLists.txt) is not set here.
    find_package(PAPI)
    if(PAPI_FOUND)
        include_directories(${PAPI_INCLUDE_DIRS})
        file(GLOB PAPI_TEST_SRCS PAPI*Tester.cpp)
    endif(PAPI_FOUND)

    # Make a list of all performance infrastructure tests we will build
    set(tests ${DUMMY_TEST_SRCS} ${COMMON_TEST_SRCS} ${PAPI_TEST_SRCS} ${OS_TEST_SRCS})

    if(CMAKE_BUILD_TYPE MATCHES "^Debug$")
        set(XOLOTL_TEST_HWCTR_DEBUGEXP 1)
        set(XOLOTL_TEST_HWCTR_EXPVAL 1150)
        set(XOLOTL_TEST_HWCTR_TOL 15)
    else(CMAKE_BUILD_TYPE MATCHES "^Debug$")
        set(XOLOTL_TEST_HWCTR_DEBUGEXP 0)
        set(XOLOTL_TEST_HWCTR_EXPVAL 1000)
        set(XOLOTL_TEST_HWCTR_TOL 5)
    endif(CMAKE_BUILD_TYPE MATCHES "^Debug$")

    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/perfTestConfig.h.in" "${CMAKE_CURRENT_BINARY_DIR}/perfTestConfig.h")

    # Define executables for each of the testers.
    foreach(test ${tests})
        message(STATUS "Defining test ${test}")
        get_filename_component(testName ${test} NAME_WE)
        add_executable(${testName} ${test})
        target_link_libraries(${testName} xolotlPerf ${MAYBE_PAPI} ${MAYBE_RT})
        add_test(${testName} ${testName})
        #add a label so the tests can be run separately
        set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})
    endforeach(test ${tests})

else(Boost_FOUND)
    message(STATUS "Not building xperf library tests; required Boost library not available")
endif(Boost_FOUND)




