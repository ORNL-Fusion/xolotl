#ifndef XCORE_FILESYSTEM_H
#define XCORE_FILESYSTEM_H

#include <xolotl/io/config.h>

#if defined(HAVE_STD_FILESYSTEM)
// We have the filesystem in the C++ standard library.
#include <filesystem>
namespace fs = std::filesystem;

#elif defined(HAVE_BOOST_FILESYSTEM)
// We have the boost filesystem.
#include "boost/filesystem.hpp"
namespace fs = boost::filesystem;

#elif defined(HAVE_STD_EXPERIMENTAL_FILESYSTEM)
// We have an experimental version of the filesystem library in the
// C++ standard library.
//
// Note: with some compilers (e.g., GCC), you may have to manually
// add a library (e.g., libstdc++fs.a for GCC) to the link line to avoid 
// undefined symbol errors when linking.
//
// TODO add checks in CMakeLists.txt to automatically find this library?
// How does the CMake cxx_check_source_compiles() command succeed
// without having found this library?
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

#else
// We have *no* filesystem library.
#error "Filesystem library (std::filesystem, std::experimental::filesystem, or Boost filesystem) required but not found."
#endif

#endif // XCORE_FILESYSTEM_H
