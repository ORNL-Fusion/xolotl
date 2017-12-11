# Set a name for this project.
PROJECT ( xolotl C CXX )
# Set the minimum required version of cmake for this project.
cmake_minimum_required(VERSION 2.8)

# Include the modules to find PETSC, ITAPS, etc. from Jed Brown
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake-modules")

# Xolotl requires C++11 to use shared pointers
INCLUDE (CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" CXX_SUPPORTS_CXX11_FLAG)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" CXX_SUPPORTS_CXX0X_FLAG)
IF (CXX_SUPPORTS_CXX11_FLAG)
    SET( CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
ELSEIF (CXX_SUPPORTS_CXX0X_FLAG)
    SET( CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")
ELSE()
    MESSAGE(FATAL_ERROR "XOLOTL requires C++11 support from C++ compiler")
ENDIF()

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

# Enable testing.
enable_testing()

# Add the preprocessor first, it is independent of the other parts of the code
# But only if Java is present
FIND_PACKAGE(Java)
IF (Java_FOUND)
    message(STATUS "The preprocessor will be built since Java was found.")
    add_subdirectory(gov.ornl.xolotl.preprocessor)
ENDIF()

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
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/diffusion)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/advection)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/temperature)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/modifiedreaction/trapmutation)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/flux)
include_directories(${CMAKE_SOURCE_DIR}/xolotlFactory/material)
include_directories(${CMAKE_SOURCE_DIR}/xolotlFactory/vizHandler)
include_directories(${CMAKE_SOURCE_DIR}/xolotlFactory/temperatureHandler)
include_directories(${CMAKE_SOURCE_DIR}/xolotlFactory/solverHandler)
include_directories(${CMAKE_SOURCE_DIR}/xolotlFactory/reactionHandler)
include_directories(${CMAKE_SOURCE_DIR}/xolotlSolver)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)
include_directories(${CMAKE_SOURCE_DIR}/xolotlViz)
include_directories(${PETSC_INCLUDES})

configure_file ("${CMAKE_SOURCE_DIR}/XolotlConfig.h.in" "${CMAKE_SOURCE_DIR}/XolotlConfig.h")

# Setup the library list
set(XOLOTL_LIBS xolotlReactants xolotlSolver xolotlIO xolotlPerf xolotlViz 
xolotlFactory xolotlCL ${PETSC_LIBRARIES} ${HDF5_LIBRARIES})

# Add an executable for the main program
ADD_EXECUTABLE (xolotl main.cpp)

# Link the reactants library
target_link_libraries(xolotl ${XOLOTL_LIBS})
