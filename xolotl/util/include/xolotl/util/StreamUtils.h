#pragma once

#include <istream>
#include <sstream>

namespace xolotl
{
namespace util
{
std::stringstream
stripComments(std::istream& is);
}
}
