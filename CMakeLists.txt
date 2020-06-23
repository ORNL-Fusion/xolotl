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

# project(xolotl C CXX)
project(xolotl CXX)

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

## Boost
# Note that we only need to list Boost component libraries that have a
# library implementation (i.e., not header only) as required components.
find_package(Boost REQUIRED COMPONENTS
    program_options
    filesystem
)

## HDF5
#set(HDF5_USE_STATIC_LIBRARIES yes)
set(HDF5_PREFER_PARALLEL ON)
enable_language(C)
find_package(HDF5 REQUIRED)
if(NOT HDF5_IS_PARALLEL)
   message(STATUS "WARNING! Detected HDF5 installation does not support parallel I/O!")
endif()

## plsm
if(NOT plsm_DIR)
    include(Getplsm)
endif()
list(APPEND CMAKE_PREFIX_PATH ${plsm_DIR})
find_package(plsm REQUIRED)

## PAPI
find_package(PAPI)

# Find VTKm
if(VTKm_DIR OR VTKM_DIR)
    # include(SetupVTKm)
    list(APPEND CMAKE_PREFIX_PATH ${VTKm_DIR} ${VTKM_DIR})
    find_package(VTKm REQUIRED)
else()
    message(STATUS "Visualization support needs explicit VTKm_DIR.")
endif()

## xolotl
add_subdirectory(xolotl)

## xconv
add_subdirectory(xconv)

include(CTest)
add_subdirectory(test)

add_subdirectory(doc)
