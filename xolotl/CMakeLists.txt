# Set the minimum required version of cmake for this project.
cmake_minimum_required(VERSION 3.4)

# Fix cmake behavior for CMAKE_CXX_STANDARD with clang on Mac OS X.
# See cmake documentation on CMP0025 policy for details.
# Must be done before doing project() command.
if(POLICY CMP0025)
    cmake_policy(SET CMP0025 NEW)
endif()

# Set policy so we can use some of the recent Jed Brown Find*.cmake files.
if(POLICY CMP0057)
    cmake_policy(SET CMP0057 NEW)
endif()

if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()

# Set a name for this project.
project(xolotl C CXX)

# Include the modules to find PETSC, ITAPS, etc. from Jed Brown
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake-modules")

# xolotl requires some C++14 features.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Tell CMake to look for static libraries
# SET(BUILD_SHARED_LIBS OFF)

find_package(MPI REQUIRED)

# Find PETSc
# We need to save and restore CMAKE_REQUIRED_{INCLUDES,LIBRARIES}
# because the FindPETSc.cmake module sets during its operation
# but does not reset them, and the leftover state then messes with
# our subsequent configuration checks.
set(SAVED_CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
set(SAVED_CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
find_package(PETSc REQUIRED)
if(PETSC_FOUND)
    include_directories(${PETSC_INCLUDES})
    set(CMAKE_REQUIRED_INCLUDES ${SAVED_CMAKE_REQUIRED_INCLUDES})
    set(CMAKE_REQUIRED_LIBRARIES ${SAVED_CMAKE_REQUIRED_LIBRARIES})
endif()

# Make sure that we only pull the static HDF5
#set(HDF5_USE_STATIC_LIBRARIES yes)
# Find HDF5 - Optional for now until we actually require it
set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 REQUIRED)
# Check for parallel I/O support in the HDF5 installation.
if(NOT HDF5_IS_PARALLEL)
   message(STATUS "WARNING! Detected HDF5 installation does not support parallel I/O!")
endif()
# Include the headers
include_directories(${HDF5_INCLUDE_DIR})

# Enable testing.
enable_testing()

# Find the Boost libraries we (potentially) use.
# Note that we only need to list Boost component libraries that have a
# library implementation (i.e., not header only) as required components.
find_package(Boost REQUIRED
    COMPONENTS
    program_options
    OPTIONAL_COMPONENTS
    chrono
    timer
    filesystem
    unit_test_framework
)

# Build the performance library first, so it can be used by other libraries.
add_subdirectory(xolotlPerf)
# Add the visualization directory
add_subdirectory(xolotlViz)
# Add the I/O, data structures and reactor directories
add_subdirectory(xolotlCore)
# Add the factory directory
add_subdirectory(xolotlFactory)
# Keep the solver for the end (it uses everything else)
add_subdirectory(xolotlSolver)

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

configure_file("${CMAKE_SOURCE_DIR}/XolotlConfig.h.in"
    "${CMAKE_BINARY_DIR}/XolotlConfig.h"
)

# Add an executable for the main program
add_executable(xolotl main.cpp)
target_link_libraries(xolotl
    xolotlReactants
    xolotlFactory
    xolotlSolver
    xolotlIO
    xolotlCL
    xolotlPerf
    xolotlViz
)
target_include_directories(xolotl PUBLIC ${CMAKE_BINARY_DIR})

if(Boost_UNIT_TEST_FRAMEWORK_FOUND)
    message(STATUS "Tests will be built since Boost unit test framework found.")
    add_subdirectory(tests)
endif()

add_subdirectory(../doc doc)
