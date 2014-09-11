# Set a name for this project.
PROJECT ( xolotl C CXX Java )
# Set the minimum required version of cmake for this project.
cmake_minimum_required(VERSION 2.8)

# Include the modules to find PETSC, ITAPS, etc. from Jed Brown
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake-modules")

# Find Java to run the preprocessor
FIND_PACKAGE(Java REQUIRED)
include(UseJava)
# Collect all of the preprocessor's dependencies
file(GLOB JAVA_DEPS "${CMAKE_SOURCE_DIR}/gov.ornl.xolotl.preprocessor" 
                    "${CMAKE_SOURCE_DIR}/gov.ornl.xolotl.preprocessor/deps" 
                    "${CMAKE_SOURCE_DIR}/gov.ornl.xolotl.preprocessor/deps/*" 
                    "${CMAKE_BINARY_DIR}/deps" "${CMAKE_BINARY_DIR}/deps/*")
                    
if (WIN32 AND NOT CYGWIN AND CMAKE_HOST_SYSTEM_NAME MATCHES "Windows")
     set(CMAKE_JAVA_INCLUDE_FLAG_SEP ";")
else ()
     set(CMAKE_JAVA_INCLUDE_FLAG_SEP ":")
endif()

# Add the preprocessor dependencies to the Java include path
foreach(dep ${JAVA_DEPS})
    set(CMAKE_JAVA_INCLUDE_PATH_FINAL "${dep}${CMAKE_JAVA_INCLUDE_FLAG_SEP}${CMAKE_JAVA_INCLUDE_PATH_FINAL}")
endforeach(dep ${JAVA_DEPS})
# Create a directory to contain the compiled preprocessor Java files
FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/preprocessor)
set(PREPROCESSOR_DIR ${CMAKE_CURRENT_BINARY_DIR}/preprocessor)
# Tell CMake where to look for the preprocessor source files
set(PREPROCESSOR_SRC ${CMAKE_SOURCE_DIR}/gov.ornl.xolotl.preprocessor/src/gov/ornl/xolotl/preprocessor)
set(CMAKE_JAVA_TARGET_OUTPUT_DIR ${PREPROCESSOR_DIR})
set(CMAKE_JAVA_CLASS_OUTPUT_PATH ${PREPROCESSOR_DIR})

# Compile the preprocessor files
add_jar(xolotlPreprocessor ${PREPROCESSOR_SRC}/Main.java ${PREPROCESSOR_SRC}/Preprocessor.java ${PREPROCESSOR_SRC}/Arguments.java 
        ${PREPROCESSOR_SRC}/Cluster.java ${PREPROCESSOR_SRC}/BindingEnergyEngine.java 
        OUTPUT_DIR ${CMAKE_JAVA_TARGET_OUTPUT_DIR})
get_target_property(_jarFile xolotlPreprocessor JAR_FILE)
get_target_property(_classDir xolotlPreprocessor CLASSDIR)
# The preprocessor .class files will be written to preprocessor/CMakeFiles/xolotlPreprocessor.dir/gov/ornl/xolotl/preprocessor
message(STATUS "xolotlPreprocessor class compiled to ${_classDir}")

# Enable testing.
enable_testing()

if ( ( DEFINED ENV{JUNIT_JAR} ) AND (DEFINED ENV{HAMCREST_JAR} ) )
    message(STATUS "Building Java JUnit tests since JUNIT_JAR and HAMCREST_JAR are defined.")

    include(JUnit)

    # Tell CMake where to look for the junit.jar file, JUNIT_JAR needs to be added to your path
    set(JUNIT_JAR $ENV{JUNIT_JAR})
    set(CMAKE_JAVA_INCLUDE_PATH_FINAL "${JUNIT_JAR}${CMAKE_JAVA_INCLUDE_FLAG_SEP}${CMAKE_JAVA_INCLUDE_PATH_FINAL}")

    # Tell CMake where to look for the hamcrest .jar file, HAMCREST_JAR needs to be added to your path
    set(HAMCREST_JAR $ENV{HAMCREST_JAR})
    set(CMAKE_JAVA_INCLUDE_PATH_FINAL "${HAMCREST_JAR}${CMAKE_JAVA_INCLUDE_FLAG_SEP}${CMAKE_JAVA_INCLUDE_PATH_FINAL}")

    # Collect all of the preprocessor's tests
    file(GLOB JAVA_JUNIT "${CMAKE_SOURCE_DIR}/gov.ornl.xolotl.preprocessor.test/src/gov/ornl/xolotl/preprocessor/test/*.java")                  

    file(GLOB JAVA_FILES "${PREPROCESSOR_SRC}/*.java")
    # Add the preprocessor dependencies to the Java include path
    foreach(jfile ${JAVA_FILES})
        set(CMAKE_JAVA_INCLUDE_PATH_FINAL "${JFILE}${CMAKE_JAVA_INCLUDE_FLAG_SEP}${CMAKE_JAVA_INCLUDE_PATH_FINAL}")
    endforeach(jfile ${JAVA_FILES})

    # Add the preprocessor tests to the Java include path
    foreach(test ${JAVA_JUNIT})
        set(CMAKE_JAVA_INCLUDE_PATH_FINAL "${CMAKE_JAVA_INCLUDE_PATH_FINAL}${CMAKE_JAVA_INCLUDE_FLAG_SEP}${test}")
    endforeach(test ${JAVA_JUNIT})

    # Create a directory to contain the compiled preprocessor Java test files
    FILE(MAKE_DIRECTORY ${PREPROCESSOR_DIR}/test)
    set(PREPROCESSOR_TEST_DIR ${PREPROCESSOR_DIR}/test)
    # Tell CMake where to look for the preprocessor's test source files
    set(PREPROCESSOR_TEST_SRC ${CMAKE_SOURCE_DIR}/gov.ornl.xolotl.preprocessor.test/src/gov/ornl/xolotl/preprocessor/test)
    set(CMAKE_JAVA_TARGET_OUTPUT_DIR ${PREPROCESSOR_TEST_DIR})
    set(CMAKE_JAVA_CLASS_OUTPUT_PATH ${PREPROCESSOR_TEST_DIR})

    set(CMAKE_JAVA_INCLUDE_PATH
       ${JUNIT_JAR}
       ${HAMCREST_JAR}
       ${CMAKE_JAVA_INCLUDE_PATH_FINAL}
    )
   
    # Compile the preprocessor tests
    message(STATUS "Making preprocessor test ${PREPROCESSOR_TEST_SRC}/BindingEnergyEngineTest.java")
    add_jar(BindingEnergyEngineTest ${PREPROCESSOR_TEST_SRC}/BindingEnergyEngineTest.java ${PREPROCESSOR_SRC}/BindingEnergyEngine.java 
        OUTPUT_DIR ${CMAKE_JAVA_TARGET_OUTPUT_DIR})
    get_target_property(_jarFile BindingEnergyEngineTest JAR_FILE)
    # The preprocessor .class files will be written to /preprocessor/test/CMakeFiles/BindingEnergyEngineTest.dir
    get_target_property(_classDir BindingEnergyEngineTest CLASSDIR)
    message(STATUS "BindingEnergyEngineTest class compiled to ${_classDir}")
    # Define the preprocessor JUnit test target
    add_junit_test(BindingEnergyEngineTest 
        CLASSPATH 
            ${CMAKE_JAVA_INCLUDE_PATH} 
            ${_classDir}
        TESTS 
            gov.ornl.xolotl.preprocessor.test.BindingEnergyEngineTest
    )


    # Compile the preprocessor tests
    message(STATUS "Making preprocessor test ${PREPROCESSOR_TEST_SRC}/PreprocessorTest.java")
    add_jar(PreprocessorTest ${PREPROCESSOR_TEST_SRC}/PreprocessorTest.java ${PREPROCESSOR_SRC}/Preprocessor.java ${PREPROCESSOR_SRC}/Arguments.java 
        ${PREPROCESSOR_SRC}/Cluster.java ${PREPROCESSOR_SRC}/BindingEnergyEngine.java 
        OUTPUT_DIR ${CMAKE_JAVA_TARGET_OUTPUT_DIR})
    get_target_property(_jarFile PreprocessorTest JAR_FILE)
    # The preprocessor test .class files will be written to preprocessor/test/CMakeFiles/PreprocessorTest.dir 
    get_target_property(_classDir PreprocessorTest CLASSDIR)
    message(STATUS "PreprocessorTest class compiled to ${_classDir}")
    # Define the preprocessor JUnit test target
    add_junit_test(PreprocessorTest 
        CLASSPATH 
            ${CMAKE_JAVA_INCLUDE_PATH} 
            ${_classDir}
        TESTS 
            gov.ornl.xolotl.preprocessor.test.PreprocessorTest
    )

    # Compile the preprocessor tests
    message(STATUS "Making preprocessor test ${PREPROCESSOR_TEST_SRC}/ArgumentsTest.java")
    add_jar(ArgumentsTest ${PREPROCESSOR_TEST_SRC}/ArgumentsTest.java ${PREPROCESSOR_SRC}/Arguments.java 
        OUTPUT_DIR ${CMAKE_JAVA_TARGET_OUTPUT_DIR})
    get_target_property(_jarFile ArgumentsTest JAR_FILE)
    # The preprocessor .class files will be written to /preprocessor/test/CMakeFiles/ArgumentsTest.dir
    get_target_property(_classDir ArgumentsTest CLASSDIR)
    message(STATUS "ArgumentsTest class compiled to ${_classDir}")
    # Define the preprocessor JUnit test target 
    add_junit_test(ArgumentsTest 
        CLASSPATH 
            ${CMAKE_JAVA_INCLUDE_PATH} 
            ${_classDir}
        TESTS 
            gov.ornl.xolotl.preprocessor.test.ArgumentsTest
    )

    # Compile the preprocessor tests
    message(STATUS "Making preprocessor test ${PREPROCESSOR_TEST_SRC}/ClusterTest.java")
    add_jar(ClusterTest ${PREPROCESSOR_TEST_SRC}/ClusterTest.java ${PREPROCESSOR_SRC}/Cluster.java 
        OUTPUT_DIR ${CMAKE_JAVA_TARGET_OUTPUT_DIR})
    get_target_property(_jarFile ClusterTest JAR_FILE)
    # The preprocessor .class files will be written to preprocessor/test/CMakeFiles/ClusterTest.dir
    get_target_property(_classDir ClusterTest CLASSDIR)
    message(STATUS "ClusterTest class compiled to ${_classDir}")
    # Define the preprocessor JUnit test target
    add_junit_test(ClusterTest 
        CLASSPATH 
            ${CMAKE_JAVA_INCLUDE_PATH} 
            ${_classDir}
        TESTS 
            gov.ornl.xolotl.preprocessor.test.ClusterTest
    )

else ()
    message(STATUS "The preprocessor JUnit tests will not be built since JUNIT_JAR and HAMCREST_JAR are not defined.")
    message(STATUS " CMake needs to know where to find junit.jar and the hamcrest jar in order to build these tests.")
    
endif ()

# Xolotl requires C++11 to use shared pointers
list( APPEND CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")

# Tell CMake to look for static libraries
SET(BUILD_SHARED_LIBS OFF)

# Find PETSc
FIND_PACKAGE(PETSc REQUIRED)

# Make sure that we only pull the static HDF5
#set(HDF5_USE_STATIC_LIBRARIES yes)
# Find HDF5 - Optional for now until we actually require it
FIND_PACKAGE(HDF5 REQUIRED)
# Check for parallel I/O support in the HDF5 installation.
IF (NOT HDF5_IS_PARALLEL)
   message(STATUS "WARNING! Detected HDF5 installation does not support parallel I/O!")
ENDIF()
# Include the headers
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})

# Build the performance library first, so it can be used by other libraries.
add_subdirectory(xolotlPerf)
# Add the visualization directory
add_subdirectory(xolotlViz)
# Add the I/O, data structures and reactor directories
add_subdirectory(xolotlCore)
add_subdirectory(xolotlSolver)

# Add the tests is Boost is available
FIND_PACKAGE(Boost)
IF (Boost_FOUND)
    message(STATUS "Boost version ${Boost_VERSION} found.")
    message(STATUS "Tests will be built since Boost was found.")
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    add_subdirectory(tests)
ENDIF()

# Report package information
message(STATUS "----- Configuration Information -----")
message(STATUS "PETSc includes = ${PETSC_INCLUDES}")
message(STATUS "PETSc libraries = ${PETSC_LIBRARIES}")
message(STATUS "----- End Configuration Information -----")

# Create the documentation directory and set up the doc target.
set(DOC_DIR "documentation")
if(NOT EXISTS ${DOC_DIR})
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${DOC_DIR}")
endif(NOT EXISTS ${DOC_DIR})
set(DOXYFILE_OUTPUT_DIR "${DOC_DIR}")
set(DOXYFILE_LATEX "ON")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
include(UseDoxygen OPTIONAL)

# Include directories for the main program
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/io)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/commandline)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters)
include_directories(${CMAKE_SOURCE_DIR}/xolotlSolver)
include_directories(${CMAKE_SOURCE_DIR}/xolotlSolver/temperature)
include_directories(${CMAKE_SOURCE_DIR}/xolotlSolver/flux)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)
include_directories(${CMAKE_SOURCE_DIR}/xolotlViz)
include_directories(${PETSC_INCLUDES})

configure_file ("${CMAKE_SOURCE_DIR}/XolotlConfig.h.in" "${CMAKE_SOURCE_DIR}/XolotlConfig.h")

# Setup the library list
set(XOLOTL_LIBS xolotlReactants xolotlSolver xolotlIO xolotlPerf xolotlViz xolotlCL ${PETSC_LIBRARIES} ${HDF5_LIBRARIES})

# Add an executable for the main program
ADD_EXECUTABLE (xolotl main.cpp)

# Link the reactants library
target_link_libraries(xolotl ${XOLOTL_LIBS})