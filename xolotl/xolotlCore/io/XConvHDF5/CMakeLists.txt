
if (HAVE_FILESYSTEM_LIBRARY AND Boost_PROGRAM_OPTIONS_FOUND)

    message(STATUS "Building XConvHDF5 utility")

    # Define the executable and the sources needed to build it.
    add_executable(xconv main.cpp)

    # Specify libraries needed to build the executable.
    target_link_libraries(xconv xolotlIO ${Boost_LIBRARIES})
else()
    message(STATUS "Not building XConvHDF5 utility (have required Boost libraries)")
endif()

