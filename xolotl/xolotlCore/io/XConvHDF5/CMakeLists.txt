
# Find the Boost libraries we (potentially) use.
# Note that we only need to list Boost component libraries that have a
# library implementation (i.e., not header only) as required components.
FIND_PACKAGE(Boost OPTIONAL_COMPONENTS program_options filesystem)
IF (Boost_FOUND)
    message(STATUS "Boost version ${Boost_VERSION} found.")
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ENDIF()
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/io/)

# Define the executable and the sources needed to build it.
add_executable(xconv main.cpp)

# Specify libraries needed to build the executable.
target_link_libraries(xconv xolotlIO ${Boost_LIBRARIES})

