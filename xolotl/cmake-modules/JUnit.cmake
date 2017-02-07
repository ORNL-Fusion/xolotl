# This CMake file was taken, then modified from https://github.com/glejeune/cmake-samples/blob/master/sample8/cmake/JUnit.cmake
#
# This CMake module provides support for JUnit testing
#
#   add_junit_test(<target name> 
#       CLASSPATH [path1 ...]
#       TESTS [class1 ...]
#   )

include(FindJava)

function(add_junit_test TARGET_NAME)

   if (WIN32 AND NOT CYGWIN)
      set(SEPARATOR ";")
   else (WIN32 AND NOT CYGWIN)
      set(SEPARATOR ":")
   endif(WIN32 AND NOT CYGWIN)

   foreach (ARG ${ARGN})
      if (ARG MATCHES "CLASSPATH" OR ARG MATCHES "TESTS")
         set(TYPE ${ARG})

      else ()

         if (TYPE MATCHES "CLASSPATH")
            set(CLASSPATH "${CLASSPATH}${SEPARATOR}${ARG}")

         elseif (TYPE MATCHES "TESTS")
            set(TESTS ${TESTS} ${ARG})

         endif()

      endif()

   endforeach(ARG)

   add_test(NAME ${TARGET_NAME} 
       COMMAND 
           ${Java_JAVA_EXECUTABLE} 
           -Djava.library.path=${CMAKE_SOURCE_DIR}/gov.ornl.xolotl.preprocessor/deps
           -classpath ${CLASSPATH} org.junit.runner.JUnitCore ${TESTS}
   )

endfunction(add_junit_test)