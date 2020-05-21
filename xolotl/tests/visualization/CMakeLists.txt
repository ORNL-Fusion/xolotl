#Set the package name
SET(PACKAGE_NAME "xolotl.tests.viz")

#Set the description
SET(PACKAGE_DESCRIPTION "Tests for the Xolotl visualization package")

# Find EAVL
#FIND_PACKAGE(EAVL)

# Find MESA
#FIND_PACKAGE(MESA)

# Find VTKm
if(VTKM_DIR)
    include(${CMAKE_SOURCE_DIR}/cmake-modules/SetupVTKm.cmake)
else()
    message(STATUS "Visualization support needs explicit VTKm_DIR.")
endif()

#Get the test files
if(VTKM_FOUND)
    file(GLOB tests *.cpp)
else(VTKM_FOUND)
    file(GLOB tests DummyPlotTester.cpp DummyDataProviderTester.cpp)
endif(VTKM_FOUND)

#Make executables and link libraries for testers
foreach(test ${tests})
    message(STATUS "Making test ${test}")
    get_filename_component(testName ${test} NAME_WE)
    add_executable(${testName} ${test})
    target_link_libraries(${testName} xolotlViz)
    add_test(${testName} ${testName})
    #add a label so the tests can be run separately
    set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})
endforeach(test ${tests})
