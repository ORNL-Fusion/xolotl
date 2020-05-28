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

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(xolotl C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")

## MPI
find_package(MPI REQUIRED)

## PETSc
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
# Report package information
message(STATUS "----- Configuration Information -----")
message(STATUS "PETSc includes = ${PETSC_INCLUDES}")
message(STATUS "PETSc libraries = ${PETSC_LIBRARIES}")
message(STATUS "----- End Configuration Information -----")

## Boost
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

# We need a filesystem library.
# We can use one of several such libraries (because the APIs are so similar).
# We would like to just check for a symbol or a header,
# but that isn't sufficient with some recent C++ compilers,
# either because they put the contents in an unexpected namespace
# (e.g., clang++ 7 putting it in a filesystem other than std::filesystem
# when compiling using the C++11 standard) or because the
# compiler can't handle the CMake tests for symbols.
include(CheckCXXSourceCompiles)
message(STATUS "Looking for C++ std::filesystem")
check_cxx_source_compiles("#include <filesystem>
int main() { std::filesystem::path test; return 0; }
" STD_FILESYSTEM_COMPILES)
if(NOT STD_FILESYSTEM_COMPILES)
    message(STATUS "Code using std::filesystem fails to compile/link.")
    set(HAVE_STD_FILESYSTEM OFF)
else()
    set(HAVE_STD_FILESYSTEM ON)
endif()

message(STATUS "Looking for C++ std::experimental::filesystem")
check_cxx_source_compiles("#include <experimental/filesystem>
int main() { std::experimental::filesystem::path test; test.has_filename(); return 0; }
" STD_EXPERIMENTAL_FILESYSTEM_COMPILES)
if(NOT STD_EXPERIMENTAL_FILESYSTEM_COMPILES)
    message(STATUS "Code using std::experimental::filesystem fails to compile/link.")
    set(HAVE_STD_EXPERIMENTAL_FILESYSTEM OFF)
else()
    set(HAVE_STD_EXPERIMENTAL_FILESYSTEM ON)
endif()

# Find the Boost libraries we (potentially) use.
# Note that we only need to list Boost component libraries that have a
# library implementation (i.e., not header only) as required components.
if (Boost_FILESYSTEM_FOUND)
    message(STATUS "Have Boost filesystem library.")
    set(HAVE_BOOST_FILESYSTEM ON)
endif()

# Set a single variable indicating whether we have a filesystem library.
if(HAVE_STD_FILESYSTEM OR HAVE_STD_EXPERIMENTAL_FILESYSTEM OR Boost_FILESYSTEM_FOUND)
    set(HAVE_FILESYSTEM_LIBRARY ON)
else()
    set(HAVE_FILESYSTEM_LIBRARY OFF)
endif()


## HDF5
#set(HDF5_USE_STATIC_LIBRARIES yes)
set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 REQUIRED)
if(NOT HDF5_IS_PARALLEL)
   message(STATUS "WARNING! Detected HDF5 installation does not support parallel I/O!")
endif()

## PAPI
find_package(PAPI)

## plsm
list(APPEND CMAKE_PREFIX_PATH ${plsm_DIR})
find_package(plsm REQUIRED)

## xolotl
add_subdirectory(xolotl)

## xconv
add_subdirectory(xconv)

include(CTest)
# add_subdirectory(test)

# add_subdirectory(doc)
