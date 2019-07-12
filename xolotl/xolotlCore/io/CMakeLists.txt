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
