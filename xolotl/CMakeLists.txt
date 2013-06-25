# Set a name for this project.
PROJECT ( xolotl C CXX )
# Set the minimum required version of cmake for this project.
cmake_minimum_required(VERSION 2.8)

# Include the modules to find PETSC, ITAPS, etc. from Jed Brown
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake-modules")

# Xolotl requires C++11 to use shared pointers
list( APPEND CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")

# Find Boost
FIND_PACKAGE(Boost REQUIRED)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

# Find PETSc
FIND_PACKAGE(PETSc REQUIRED)

# Enable testing.
enable_testing()

# Add the I/O, data structures and reactor directories
add_subdirectory(xolotlCore)
add_subdirectory(xolotlSolver)
add_subdirectory(tests)

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
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters)
include_directories(${CMAKE_SOURCE_DIR}/xolotlSolver)
include_directories(${PETSC_INCLUDES})

# Setup the library list
set(XOLOTL_LIBS xolotlReactants xolotlSolver xolotlIO ${PETSC_LIBRARIES})

# Add an executable for the main program
ADD_EXECUTABLE ( xolotl main.cpp)

# Link the reactants library
target_link_libraries(xolotl ${XOLOTL_LIBS})