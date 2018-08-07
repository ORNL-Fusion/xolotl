#ifndef XCORE_FILESYSTEM_H
#define XCORE_FILESYSTEM_H

#include "xolotlCore/io/XDMConfig.h"

#if defined(HAVE_STD_FILESYSTEM)
// We have the filesystem in the C++ standard library.
#include <filesystem>
namespace fs = std::filesystem;

#elif defined(HAVE_STD_EXPERIMENTAL_FILESYSTEM)
// We have an experimental version of the filesystem library in the
// C++ standard library.
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

#elif defined(HAVE_BOOST_FILESYSTEM)
// We have the boost filesystem.
#include "boost/filesystem.hpp"
namespace fs = boost::filesystem;

#else
// We have *no* filesystem library.
#error "Filesystem library (std::filesystem, std::experimental::filesystem, or Boost filesystem) required but not found."
#endif

#endif // XCORE_FILESYSTEM_H
