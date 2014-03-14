#Set the package name
SET(PACKAGE_NAME "xolotl.tests.performance")

#Set the description
SET(PACKAGE_DESCRIPTION "Tests for the Xolotl performance package")

# Find GPTL
FIND_PACKAGE(GPTL)

# Find PAPI
FIND_PACKAGE(PAPI)

if(GPTL_FOUND AND PAPI_FOUND)

#Include directories from the source and boost binaries
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/dummy)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/standard)
include_directories(${Boost_BINARY_DIRS}) 
include_directories(${GPTL_INCLUDE_DIRS})
include_directories(${PAPI_INCLUDE_DIRS})

#Link the directories
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlPerf)
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlPerf/dummy)
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlPerf/standard)
    
    #Get the test files
    file(GLOB tests *Tester.cpp)
    
#If boost was found, create tests
if(Boost_FOUND)
    #Make executables and link libraries for testers
    foreach(test ${tests})
        message(STATUS "Making test ${test}")
        get_filename_component(testName ${test} NAME_WE)
        add_executable(${testName} ${test})
        target_link_libraries(${testName} xolotlPerf ${GPTL_LIBRARIES} ${PAPI_LIBRARIES})
        add_test(${testName} ${testName}) 
        #add a label so the tests can be run separately
        set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})   
    endforeach(test ${tests})
endif(Boost_FOUND) 

else(GPTL_FOUND AND PAPI_FOUND)

#Include directories from the source and boost binaries
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/dummy)
include_directories(${Boost_BINARY_DIRS}) 

#Link the directories
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlPerf)
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlPerf/dummy)
    
    #Get the test files
 #   file(GLOB tests *Tester.cpp dummy/*Tester.cpp standard/*Tester.cpp) 
    file(GLOB tests DummyTimerTester.cpp DummyEventCounterTester.cpp DummyHardwareCounterTester.cpp)
    
#If boost was found, create tests
if(Boost_FOUND)
    #Make executables and link libraries for testers
    foreach(test ${tests})
        message(STATUS "Making test ${test}")
        get_filename_component(testName ${test} NAME_WE)
        add_executable(${testName} ${test})
        target_link_libraries(${testName} xolotlPerf)
        add_test(${testName} ${testName}) 
        #add a label so the tests can be run separately
        set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})   
    endforeach(test ${tests})
    endif(Boost_FOUND)

endif(GPTL_FOUND AND PAPI_FOUND)





