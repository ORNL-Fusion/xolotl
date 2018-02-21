# This CMake file takes care of compiling the Java preprocessor

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
        ${PREPROCESSOR_SRC}/Cluster.java ${PREPROCESSOR_SRC}/FormationEnergyEngine.java 
        OUTPUT_DIR ${CMAKE_JAVA_TARGET_OUTPUT_DIR})
get_target_property(_jarFile xolotlPreprocessor JAR_FILE)
get_target_property(_classDir xolotlPreprocessor CLASSDIR)
# The preprocessor .class files will be written to preprocessor/CMakeFiles/xolotlPreprocessor.dir/gov/ornl/xolotl/preprocessor
message(STATUS "xolotlPreprocessor class compiled to ${_classDir}")

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
   
    # Compile the formation energy engine tests
    message(STATUS "Making preprocessor test ${PREPROCESSOR_TEST_SRC}/FormationEnergyEngineTest.java")
    add_jar(FormationEnergyEngineTest ${PREPROCESSOR_TEST_SRC}/FormationEnergyEngineTest.java ${PREPROCESSOR_SRC}/FormationEnergyEngine.java 
        OUTPUT_DIR ${CMAKE_JAVA_TARGET_OUTPUT_DIR})
    get_target_property(_jarFile FormationEnergyEngineTest JAR_FILE)
    # The preprocessor .class files will be written to /preprocessor/test/CMakeFiles/FormationEnergyEngineTest.dir
    get_target_property(_classDir FormationEnergyEngineTest CLASSDIR)
    message(STATUS "FormationEnergyEngineTest class compiled to ${_classDir}")
    # Define the preprocessor JUnit test target
    add_junit_test(FormationEnergyEngineTest 
        CLASSPATH 
            ${CMAKE_JAVA_INCLUDE_PATH} 
            ${_classDir}
        TESTS 
            gov.ornl.xolotl.preprocessor.test.FormationEnergyEngineTest
    )


    # Compile the preprocessor tests
    message(STATUS "Making preprocessor test ${PREPROCESSOR_TEST_SRC}/PreprocessorTest.java")
    add_jar(PreprocessorTest ${PREPROCESSOR_TEST_SRC}/PreprocessorTest.java ${PREPROCESSOR_SRC}/Preprocessor.java ${PREPROCESSOR_SRC}/Arguments.java 
        ${PREPROCESSOR_SRC}/Cluster.java ${PREPROCESSOR_SRC}/FormationEnergyEngine.java 
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

    # Compile the arguments tests
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

else ()
    message(STATUS "The preprocessor JUnit tests will not be built since JUNIT_JAR and HAMCREST_JAR are not defined.")
    message(STATUS "CMake needs to know where to find junit.jar and the hamcrest jar in order to build these tests.")
    
endif ()