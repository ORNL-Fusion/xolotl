if(HAVE_FILESYSTEM_LIBRARY AND Boost_PROGRAM_OPTIONS_FOUND)
    message(STATUS "Building XConvHDF5 utility")
    add_executable(xconv main.cpp)
    target_link_libraries(xconv xolotlIO Boost::program_options)
else()
    message(STATUS "Not building XConvHDF5 utility (don't have required Boost libraries)")
endif()

