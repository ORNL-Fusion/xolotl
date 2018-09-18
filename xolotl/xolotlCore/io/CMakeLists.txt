#Set the package name
SET(PACKAGE_NAME "xolotl.io")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl I/O")
#Set the library name
SET(LIBRARY_NAME "xolotlIO")

# Define our library and specify the sources needed to build it.
add_library(${LIBRARY_NAME} STATIC
            HDF5File.cpp
            HDF5FileType.cpp
            HDF5FileAttribute.cpp
            HDF5FileGroup.cpp
            HDF5FileDataSpace.cpp
            HDF5FileDataSet.cpp
            XFile.cpp
            MPIUtils.cpp)

# We need a filesystem library.
# We can use one of several (because the APIs are so similar).
# Note: in the following, we'd like to use CheckCXXSymbolExists(),
# but as of CMake-3.10 it seems to do its test in a way that
# is not accepted by some recent C++ compilers (e.g., clang++ 5.0).
# So we use try_compile() instead.
include(CheckIncludeFileCXX)
check_include_file_cxx("filesystem" HAVE_STD_FILESYSTEM)
check_include_file_cxx("experimental/filesystem" HAVE_STD_EXPERIMENTAL_FILESYSTEM)
if(HAVE_STD_EXPERIMENTAL_FILESYSTEM)
    # In some environments (e.g., OS X), the libc++ that is found
    # by default does not have the experimental filesystem support
    # built into it.  Instead, we have to link explicitly against
    # libc++experimental.  See whether we can use *that* library,
    # and if not, alter our notion of whether std::experimental::filesystem
    # is usable.
    find_library(HAVE_STD_EXPERIMENTAL_FILESYSTEM_LIBRARY c++experimental)
    if(NOT HAVE_STD_EXPERIMENTAL_FILESYSTEM_LIBRARY)
        message(STATUS "Can't find the library implementing std::experimental::filesystem.")
        set(HAVE_STD_EXPERIMENTAL_FILESYSTEM OFF)
    endif()
endif()

# Find the Boost libraries we (potentially) use.
# Note that we only need to list Boost component libraries that have a
# library implementation (i.e., not header only) as required components.
if(Boost_FILESYSTEM_LIBRARY_RELEASE)
    message(STATUS "Already found Boost filesystem library.")
    set(HAVE_BOOST_FILESYSTEM ON)
    include_directories(${Boost_INCLUDE_DIR})
endif()

#Collect all header filenames in this project 
#and glob them in HEADERS
# TODO we don't need to install anything when building this internal library,
# so no need to figure out headers.
#file(GLOB HEADERS *.h)

# Where to find include files.
include_directories(${CMAKE_SOURCE_DIR}
                    ${CMAKE_BINARY_DIR}
                    ${MPI_INCLUDE_PATH}
                    ${HDF5_INCLUDE_DIR}
                    ${CMAKE_SOURCE_DIR}/xolotlCore/
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters/
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/feclusters/
                    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants/neclusters/)
target_link_libraries(${LIBRARY_NAME} ${MPI_LIBRARIES} ${HDF5_LIBRARIES} ${Boost_LIBRARIES})

add_subdirectory(XConvHDF5)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION ${PROJECT_BINARY_DIR}/include)

# Configure file for the I/O library.
# Note: this must be done after all other checks, or the contents of the
# generated file doesn't take the check results into account.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/XDMConfig.h.in" "${CMAKE_CURRENT_BINARY_DIR}/XDMConfig.h")
