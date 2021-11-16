cmake_minimum_required(VERSION 3.13)

# Fix cmake behavior for CMAKE_CXX_STANDARD with clang on Mac OS X.
# See cmake documentation on CMP0025 policy for details.
# Must be done before doing project() command.
if(POLICY CMP0025)
    cmake_policy(SET CMP0025 NEW)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
include(VersionFromGit)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(Xolotl
    VERSION ${VERSION}
    LANGUAGES CXX
)

## MPI
find_package(MPI REQUIRED)

## PETSc
find_package(PETSc REQUIRED)
option(Xolotl_USE_64BIT_INDEX_TYPE "" ${PETSc_USE_64BIT_INDICES})

## Boost
# Note that we only need to list Boost component libraries that have a
# library implementation (i.e., not header only) as required components.
find_package(Boost REQUIRED COMPONENTS
    filesystem
    log_setup
    log
    program_options
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
include(SetupPlsm)
list(APPEND CMAKE_PREFIX_PATH ${plsm_DIR})
find_package(plsm REQUIRED)

## PAPI
find_package(PAPI)

## VTKm
if(VTKm_DIR OR VTKM_DIR)
    list(APPEND CMAKE_PREFIX_PATH ${VTKm_DIR} ${VTKM_DIR})
    find_package(VTKm REQUIRED)
else()
    message(STATUS "Visualization support needs explicit VTKm_DIR.")
endif()

include(CTest)
include(XolotlExtraTesting)

## xolotl
set(XOLOTL_USE_64BIT_INDEX_TYPE ${Xolotl_USE_64BIT_INDEX_TYPE} CACHE INTERNAL "")
add_subdirectory(xolotl)

## xconv
add_subdirectory(xconv)

add_subdirectory(test)

add_subdirectory(doc)
