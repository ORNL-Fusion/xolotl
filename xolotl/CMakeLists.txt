#Set a name for this project.
PROJECT ( xolotl C CXX )
#Set the minimum required version of cmake for this project.
cmake_minimum_required(VERSION 2.8)

#Xolotl requires C++11 to use shared pointers
list( APPEND CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")

# Add an executabe
ADD_EXECUTABLE ( xolotl main.cpp)

#Find boost customarily
FIND_PACKAGE(Boost REQUIRED)

#Enable testing.
enable_testing()

#Add the I/O, data structures and reactor directories
add_subdirectory(xolotlCore)
add_subdirectory(tests)

#Create the documentation directory and set up the doc target.
set(DOC_DIR "documentation")
if(NOT EXISTS ${DOC_DIR})
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "${DOC_DIR}")
endif(NOT EXISTS ${DOC_DIR})
set(DOXYFILE_OUTPUT_DIR "${DOC_DIR}")
set(DOXYFILE_LATEX "ON")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
include(UseDoxygen OPTIONAL)

